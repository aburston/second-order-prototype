"""Map the three level prototype against Van der Pol from mu = 0.1 to 5.

Run ``python3 campaign.py all`` (hours) or one stage at a time::

    python3 campaign.py survey     # Van der Pol's plateau edges at each mu
    python3 campaign.py fit        # fit the model at each mu, in priority order
    python3 campaign.py verify     # sweep each fitted model beside Van der Pol
    python3 campaign.py formula    # power laws through the fitted parameters

Every result is appended to ``campaign/results.json`` as it completes and
committed, so the campaign can be stopped and resumed and the document
written from partial data. ``THREELEVEL.md`` reads its tables from that
file through ``figures.fig_campaign`` and ``campaign.report``.

The recipe is the one ``THREELEVEL.md`` describes: at drive amplitude 5,
fit the model's three ratios and two edges so that the end of its 1:1
plateau and the start and end of its 3:1 plateau sit on Van der Pol's,
with the free cycle held within 20%. The fits run from the middle of the
range outward, each started from the interpolation of the fits already
made and, once four exist, from the power laws through them.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

import staircase
import vanderpol

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "campaign", "results.json")
MUS = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
#: Fit order: between the known points first, then outward.
PRIORITY = (2.0, 3.0, 0.5, 1.5, 4.0, 0.3, 0.7, 0.2, 0.1, 1.0, 5.0)
KNOWN = {1.0: staircase.THREE_FITTED_MU1, 5.0: staircase.THREE_FITTED}
AMP = staircase.CMP_AMP


def load():
    if os.path.exists(RESULTS):
        return json.load(open(RESULTS))
    return {"targets": {}, "fits": {}, "verify": {}, "formula": {}}


def save(res, message):
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    json.dump(res, open(RESULTS, "w"), indent=1, sort_keys=True)
    try:
        subprocess.run(["git", "add", RESULTS], cwd=HERE, check=True)
        subprocess.run(["git", "commit", "-q", "-m", message + "\n\n"
                        "Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>\n"
                        "Claude-Session: https://claude.ai/code/session_01TCkveJPUbZ33KP7WqV2MuG"],
                       cwd=HERE, check=True)
        subprocess.run(["git", "push", "-q", "origin", "HEAD"], cwd=HERE, check=True)
    except subprocess.CalledProcessError as e:
        print("  (git step failed: %s)" % e, flush=True)


def log(msg):
    print(time.strftime("%H:%M:%S "), msg, flush=True)


# ------------------------------------------------------------------ survey
def survey(res, mus=MUS):
    for mu in mus:
        key = "%g" % mu
        if key in res["targets"]:
            continue
        t0 = time.time()
        tg, labels, wl = staircase.campaign_targets(mu, amp=AMP)
        T = vanderpol.cycle(mu)[0]
        R = vanderpol.amplitude(mu)
        clean = lambda v: None if v is None else [float(v[0]), float(v[1])]
        res["targets"][key] = dict(mu=mu, w_lc=float(wl), R=float(R), T=float(T),
                                   lock1=clean(tg["lock1"]), lock3=clean(tg["lock3"]),
                                   labels=labels)
        log("survey mu=%g: R %.4f T %.3f  lock1 %s  lock3 %s  (%.0fs)"
            % (mu, R, T, tg["lock1"], tg["lock3"], time.time() - t0))
        save(res, "Campaign survey: Van der Pol's plateau edges at mu = %g" % mu)


# ------------------------------------------------------------------- fits
def power_laws(res):
    """Fit ``zeta_k = c_k mu^p_k`` and constant edges through the fits."""
    pts = sorted((float(k), v) for k, v in res["fits"].items())
    if len(pts) < 3:
        return None
    mus = np.array([m for m, _ in pts])
    lv = np.array([v["levels"] for _, v in pts])
    ed = np.array([v["edges"] for _, v in pts])
    laws = []
    for k in range(3):
        y = lv[:, k]
        sgn = np.sign(y[0])
        p, c = np.polyfit(np.log(mus), np.log(np.abs(y)), 1)
        laws.append(dict(sign=float(sgn), c=float(np.exp(c)), p=float(p)))
    return dict(levels=laws, edges=[float(e) for e in ed.mean(axis=0)],
                edge_spread=[float(e) for e in ed.std(axis=0)], n=len(pts))


def predict(law, mu):
    lv = tuple(l["sign"]*l["c"]*mu**l["p"] for l in law["levels"])
    return lv, tuple(law["edges"])


def start_for(res, mu):
    """Starting model: the power laws if four fits exist, else interpolation
    in log mu between the nearest fitted or known models."""
    law = power_laws(res)
    if law is not None and law["n"] >= 4:
        return predict(law, mu), "power laws through %d fits" % law["n"]
    known = dict(KNOWN)
    for k, v in res["fits"].items():
        known[float(k)] = (tuple(v["levels"]), tuple(v["edges"]))
    ms = sorted(known)
    lo = max([m for m in ms if m <= mu], default=None)
    hi = min([m for m in ms if m >= mu], default=None)
    if lo is None or hi is None or lo == hi:
        a, b = (ms[0], ms[1]) if mu < ms[0] else (ms[-2], ms[-1])
    else:
        a, b = lo, hi
    (la, ea), (lb, eb) = known[a], known[b]
    w = (np.log(mu) - np.log(a))/(np.log(b) - np.log(a))
    lv = tuple(float(np.sign(x)*np.exp((1 - w)*np.log(abs(x)) + w*np.log(abs(y))))
               for x, y in zip(la, lb))
    ed = tuple(float((1 - w)*x + w*y) for x, y in zip(ea, eb))
    return (lv, ed), "log interpolation between mu = %g and %g" % (a, b)


def fit(res, mus=PRIORITY, maxfev=30):
    for mu in mus:
        key = "%g" % mu
        if key in res["fits"]:
            continue
        if key not in res["targets"]:
            survey(res, (mu,))
        tg = res["targets"][key]
        targets = {"lock1": tg["lock1"], "lock3": tg["lock3"]}
        if targets["lock3"] is None and targets["lock1"] is None:
            log("fit mu=%g: no plateau to fit at A = %g, skipped" % (mu, AMP))
            res["fits"][key] = dict(mu=mu, skipped="no plateau")
            save(res, "Campaign: no plateau to fit at mu = %g" % mu)
            continue
        extra = []
        if targets["lock3"] is None:
            # no 3:1 plateau at A = 5: add the 1:1 tongue's two edges at A = 1
            if "lock1_A1" not in tg:
                xt, _, _ = staircase.campaign_targets(mu, amp=1.0, r_hi=2.5)
                tg["lock1_A1"] = None if xt["lock1"] is None else [float(v) for v in xt["lock1"]]
                log("fit mu=%g: 1:1 tongue at A = 1 is %s" % (mu, tg["lock1_A1"]))
            if tg["lock1_A1"]:
                extra.append((1.0, {"lock1": tuple(tg["lock1_A1"])}))
        (lv0, ed0), how = start_for(res, mu)
        log("fit mu=%g from %s: levels %s edges %s" % (mu, how, lv0, ed0))
        t0 = time.time()
        lv, ed, found, r, T, n = staircase.fit_plateaus(
            mu, lv0, ed0, targets, tg["w_lc"], maxfev=maxfev, amp=AMP,
            log=lambda m: print(m, flush=True), free=(tg["R"], tg["T"]),
            extra=extra)
        res["fits"][key] = dict(mu=mu, levels=[float(z) for z in lv],
                                edges=[float(e) for e in ed],
                                found={"%s%d" % (k[0], k[1]): float(v) for k, v in found.items()},
                                extra_targets=[(a, {k: list(v) for k, v in d.items()}) for a, d in extra],
                                r=float(r), T=float(T), n_eval=int(n),
                                start=[float(v) for v in list(lv0) + list(ed0)],
                                start_how=how, seconds=time.time() - t0)
        log("fit mu=%g done: levels %s edges %s r %.4f T %.3f  %s  (%d evals, %.0fs)"
            % (mu, tuple(round(z, 3) for z in lv), tuple(round(e, 3) for e in ed),
               r, T, found, n, time.time() - t0))
        law = power_laws(res)
        if law:
            res["formula"] = law
            log("power laws: " + "  ".join("zeta%d = %+.3f mu^%.3f" % (k, l["sign"]*l["c"], l["p"])
                                            for k, l in enumerate(law["levels"]))
                + "  edges %s" % tuple(round(e, 3) for e in law["edges"]))
        save(res, "Campaign fit at mu = %g" % mu)


# ----------------------------------------------------------------- verify
def verify(res, mus=PRIORITY, r_lo=0.5, r_hi=6.0, step=0.02):
    import section
    for mu in mus:
        key = "%g" % mu
        if key in res["verify"] or key not in res["fits"] or "levels" not in res["fits"][key]:
            continue
        model = (tuple(res["fits"][key]["levels"]), tuple(res["fits"][key]["edges"]))
        wl = res["targets"][key]["w_lc"]
        oms = tuple(np.round(np.arange(r_lo, r_hi + 1e-9, step)*wl, 6))
        t0 = time.time()
        scan = staircase.system_scan(oms, {"vdp": None, "three": model}, mu=mu, amp=AMP)
        out = {}
        for tag in ("vdp", "three"):
            rows = []
            for om, lab, lam in scan[tag]:
                if lab == "chaos":
                    flow = (vanderpol.field(mu, AMP, om) if tag == "vdp"
                            else staircase.field(model[0], model[1], AMP, om))
                    y0 = list(vanderpol.cycle(mu)[1]) if tag == "vdp" else [2.0, 0.0]
                    ok, _ = section.confirm_chaos(flow, 2*np.pi/om, y0, staircase.CMP_NSKIP)
                    lab = "chaos" if ok else "torus"
                rows.append((round(om/wl, 3), lab, lam))
            out[tag] = dict(runs=[(l, float(lo), float(hi), int(n)) for l, lo, hi, n in staircase.runs(rows)],
                            chaotic=[float(r) for r, l, _ in rows if l == "chaos"])
        agree = staircase.window_agreement({"vdp": scan["vdp"], "three": scan["three"]})
        res["verify"][key] = dict(mu=mu, step=step, r_lo=r_lo, r_hi=r_hi,
                                  vdp=out["vdp"], three=out["three"],
                                  jaccard=float(agree["three"][2]), seconds=time.time() - t0)
        log("verify mu=%g: vdp chaotic %d, three chaotic %d, jaccard %.3f (%.0fs)"
            % (mu, len(out["vdp"]["chaotic"]), len(out["three"]["chaotic"]),
               agree["three"][2], time.time() - t0))
        save(res, "Campaign verification at mu = %g" % mu)


def report(res):
    print("%6s %8s %8s %8s %6s %6s %8s %8s %10s %10s %10s %10s" % (
        "mu", "zeta0", "zeta1", "zeta2", "a", "b", "r", "T", "lock1 st", "lock1 end", "lock3 st", "lock3 end"))
    for k in sorted(res["fits"], key=float):
        f = res["fits"][k]
        if "levels" not in f:
            continue
        fd = f["found"]
        print("%6s %8.3f %8.3f %8.3f %6.3f %6.3f %8.4f %8.3f %10s %10s %10s %10s" % (
            k, *f["levels"], *f["edges"], f["r"], f["T"],
            "%.3f" % fd["lock10"] if "lock10" in fd else "-",
            "%.3f" % fd["lock11"] if "lock11" in fd else "-",
            "%.3f" % fd["lock30"] if "lock30" in fd else "-",
            "%.3f" % fd["lock31"] if "lock31" in fd else "-"))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = load()
    if what in ("survey", "all"):
        survey(res)
    if what in ("fit", "all"):
        fit(res)
    if what in ("verify", "all"):
        verify(res)
    if what == "formula":
        law = power_laws(res)
        print(law)
    if what in ("report", "all"):
        report(res)
