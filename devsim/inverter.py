"""
Mixed-mode CMOS inverter: two 2D drift-diffusion devices coupled through a
circuit netlist, solved simultaneously (the DEVSIM equivalent of a
Sentaurus SDevice `System { }` mixed-mode simulation).

    VDD ─┬─────────────┐
         │          pMOS (source, body = vdd)
    in ──┼── gates      ├── out ── CL ── gnd
         │          nMOS (source, body = gnd)
    gnd ─┴─────────────┘

Analyses:
  dc    : voltage transfer characteristic (VTC), supply current
  tran  : pulse response with load capacitance -> tpHL, tpLH, tr, tf

Run:  python inverter.py dc   [--wp 1.5]
      python inverter.py tran [--wp 1.5] [--cl 10e-15]
"""

import argparse
import json
import os

import numpy as np

os.environ.setdefault("DEVSIM_MATH_LIBS", "/usr/local/lib/libmkl_rt.so.3")
import devsim as ds  # noqa: E402

import tcad  # noqa: E402

VDD = 1.8
OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def solve():
    # Circuit node voltages near 0 V make the relative error noisy (~1e-7)
    # even when absolute error is ~1e-17 V; 1e-6 relative is ~1 uV here.
    tcad.solve(rel_err=1e-6)


def build(wn_um=1.0, wp_um=1.5, cl=0.0):
    """Create both devices wired into the inverter netlist.

    width_factor = W[um] * 1e-4 converts DEVSIM's per-cm 2D currents and
    charges to a device of W um, so every circuit quantity is in real units.
    """
    for d in ds.get_device_list():
        ds.delete_device(device=d)
    for m in ds.get_mesh_list():
        ds.delete_mesh(mesh=m)
    try:
        ds.delete_circuit()
    except ds.error:
        pass  # no circuit yet

    # Netlist (same topology as the Sentaurus System{} block):
    ds.circuit_element(name="Vdd", n1="vdd", n2="0", value=0.0)
    ds.circuit_element(name="Vgnd", n1="gnd", n2="0", value=0.0)
    ds.circuit_element(name="Vin", n1="in", n2="0", value=0.0)
    ds.add_circuit_node(name="out", value=0.0, variable_update="default")
    # 100 TOhm bleed: before drift-diffusion is enabled no current flows into
    # "out", so its KCL row would be empty (singular). 18 fA at 1.8 V is
    # ~500x below the transistors' off-current, i.e. negligible.
    ds.circuit_element(name="Rbleed", n1="out", n2="0", value=1e14)
    if cl > 0:
        ds.circuit_element(name="CL", n1="out", n2="0", value=cl)

    nmap = dict(source="gnd", body="gnd", gate="in", drain="out")
    pmap = dict(source="vdd", body="vdd", gate="in", drain="out")
    pn = tcad.MOSParams("n", width_factor=wn_um * 1e-4)
    pp = tcad.MOSParams("p", width_factor=wp_um * 1e-4)
    tcad.create_device("nmos", pn, circuit_nodes=nmap)
    tcad.create_device("pmos", pp, circuit_nodes=pmap)
    solve()  # equilibrium (all sources 0)
    tcad.enable_dd("nmos", pn, circuit_nodes=nmap)
    tcad.enable_dd("pmos", pp, circuit_nodes=pmap)
    solve()


def node(n):
    return ds.get_circuit_node_value(node=n, solution="dcop")


def _apply(name, value):
    ds.circuit_alter(name=name, value=value)


def set_source(name, target, step=0.05, min_step=1e-4):
    """Ramp a source (Vdd circuit source or the Vin gate drive) adaptively."""
    cur = _state[name]
    while abs(target - cur) > 1e-12:
        dv = np.sign(target - cur) * min(step, abs(target - cur))
        try:
            _apply(name, cur + dv)
            solve()
            cur += dv
            step = min(step * 1.5, 0.1)
        except ds.error:
            _apply(name, cur)
            step /= 2
            if step < min_step:
                raise RuntimeError(f"{name} ramp failed at {cur}")
    _state[name] = cur


_state = {"Vdd": 0.0, "Vin": 0.0}


def supply_current():
    """Current delivered by VDD (A); DEVSIM reports source branch current."""
    try:
        return -ds.get_circuit_node_value(node="Vdd.I", solution="dcop")
    except ds.error:
        return float("nan")


def vtc(wn=1.0, wp=1.5, dv=0.02):
    _state.update(Vdd=0.0, Vin=0.0)
    build(wn, wp)
    set_source("Vdd", VDD)
    vin, vout, idd = [], [], []
    v = 0.0
    while v <= VDD + 1e-9:
        set_source("Vin", v, step=dv)
        vin.append(v)
        vout.append(node("out"))
        idd.append(supply_current())
        v = round(v + dv, 10)
    return np.array(vin), np.array(vout), np.array(idd)


def vtc_metrics(vin, vout):
    """Switching threshold, unity-gain points, noise margins, peak gain."""
    gain = np.gradient(vout, vin)
    vm = float(np.interp(0.0, (vout - vin)[::-1], vin[::-1]))  # Vout = Vin
    idx = np.where(gain <= -1.0)[0]
    vil = float(np.interp(-1.0, gain[: idx[0] + 1][::-1], vin[: idx[0] + 1][::-1])) if len(idx) else float("nan")
    vih = float(np.interp(-1.0, gain[idx[-1]:], vin[idx[-1]:])) if len(idx) else float("nan")
    voh = float(np.interp(vil, vin, vout))
    vol = float(np.interp(vih, vin, vout))
    return dict(
        VM=vm,
        VIL=vil,
        VIH=vih,
        VOH=voh,
        VOL=vol,
        NML=vil - vol,
        NMH=voh - vih,
        peak_gain=float(-gain.min()),
        Vout_at_Vin0=float(vout[0]),
        Vout_at_VinVDD=float(vout[-1]),
    )


def transient(wn=1.0, wp=1.5, cl=10e-15, t_rise=20e-12, t_high=300e-12, dt_fine=1e-12, dt_coarse=5e-12):
    """Input pulse 0 -> VDD -> 0; returns t, vin, vout, metrics."""
    _state.update(Vdd=0.0, Vin=0.0)
    build(wn, wp, cl)
    set_source("Vdd", VDD)
    t_total = 2 * t_rise + 2 * t_high

    def vin_at(t):
        if t < t_high / 2:
            return 0.0
        t -= t_high / 2
        if t < t_rise:
            return VDD * t / t_rise
        t -= t_rise
        if t < t_high:
            return VDD
        t -= t_high
        if t < t_rise:
            return VDD * (1 - t / t_rise)
        return 0.0

    ds.solve(type="transient_dc", absolute_error=1e10, relative_error=1e-6, maximum_iterations=40)
    edges = (t_high / 2, t_high / 2 + t_rise + t_high)  # input edge start times

    def max_step(t):
        # 1 ps from just before each input edge until the output has settled
        # (~150 ps later); 5 ps in the quiet stretches.
        for e in edges:
            if e - 10e-12 <= t <= e + 150e-12:
                return dt_fine
        return dt_coarse

    ts, vi, vo = [0.0], [0.0], [node("out")]
    t, h, n = 0.0, dt_fine, 0
    while t < t_total - 1e-18:
        h = min(h, max_step(t), t_total - t)
        _apply("Vin", vin_at(t + h))
        try:
            # charge_error=1: rely on the explicit step schedule instead of
            # DEVSIM's charge-error test (it rejected steps at the edges and
            # halved h without bound).
            ds.solve(type="transient_bdf1", absolute_error=1e10, relative_error=1e-6,
                     maximum_iterations=30, tdelta=h, charge_error=1.0)
        except ds.error:
            h /= 2
            if h < 1e-15:
                raise
            continue
        t += h
        n += 1
        ts.append(t)
        vi.append(vin_at(t))
        vo.append(node("out"))
        if n % 20 == 0:
            print(f"PROGRESS t={t*1e12:.1f} ps  Vin={vi[-1]:.3f}  Vout={vo[-1]:.4f}", flush=True)
        h = min(h * 1.5, max_step(t))
    ts, vi, vo = map(np.array, (ts, vi, vo))
    return ts, vi, vo


def crossing(t, v, level, rising, after=0.0):
    for i in range(1, len(t)):
        if t[i] <= after:
            continue
        a, b = v[i - 1], v[i]
        if (rising and a < level <= b) or (not rising and a > level >= b):
            return float(t[i - 1] + (level - a) / (b - a) * (t[i] - t[i - 1]))
    return float("nan")


def tran_metrics(t, vi, vo):
    half = VDD / 2
    tin_r = crossing(t, vi, half, True)
    tout_f = crossing(t, vo, half, False, after=tin_r - 1e-12)
    tin_f = crossing(t, vi, half, False, after=tin_r)
    tout_r = crossing(t, vo, half, True, after=tin_f - 1e-12)
    # 10-90 % transition times of the output
    t90f = crossing(t, vo, 0.9 * VDD, False, after=tin_r - 1e-12)
    t10f = crossing(t, vo, 0.1 * VDD, False, after=tin_r - 1e-12)
    t10r = crossing(t, vo, 0.1 * VDD, True, after=tin_f - 1e-12)
    t90r = crossing(t, vo, 0.9 * VDD, True, after=tin_f - 1e-12)
    tphl, tplh = tout_f - tin_r, tout_r - tin_f
    return dict(
        tpHL_ps=tphl * 1e12,
        tpLH_ps=tplh * 1e12,
        tp_avg_ps=(tphl + tplh) / 2 * 1e12,
        tfall_10_90_ps=(t10f - t90f) * 1e12,
        trise_10_90_ps=(t90r - t10r) * 1e12,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("analysis", choices=["dc", "tran", "sizing"])
    ap.add_argument("--wn", type=float, default=1.0)
    ap.add_argument("--wp", type=float, default=1.5)
    ap.add_argument("--cl", type=float, default=10e-15)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    tag = f"wp{a.wp:g}"

    if a.analysis == "dc":
        vin, vout, idd = vtc(a.wn, a.wp)
        np.savetxt(os.path.join(OUT, f"inverter_vtc_{tag}.csv"), np.c_[vin, vout, idd], delimiter=",",
                   header="Vin,Vout,Idd_A", comments="")
        m = vtc_metrics(vin, vout)
        m["peak_Idd_uA"] = float(np.nanmax(np.abs(idd)) * 1e6)
        m["Wn_um"], m["Wp_um"] = a.wn, a.wp
        with open(os.path.join(OUT, f"inverter_dc_{tag}.json"), "w") as f:
            json.dump(m, f, indent=2)
        print(json.dumps(m, indent=2))
    elif a.analysis == "tran":
        t, vi, vo = transient(a.wn, a.wp, a.cl)
        np.savetxt(os.path.join(OUT, f"inverter_tran_{tag}.csv"), np.c_[t, vi, vo], delimiter=",",
                   header="t_s,Vin,Vout", comments="")
        m = tran_metrics(t, vi, vo)
        m.update(Wn_um=a.wn, Wp_um=a.wp, CL_fF=a.cl * 1e15)
        with open(os.path.join(OUT, f"inverter_tran_{tag}.json"), "w") as f:
            json.dump(m, f, indent=2)
        print(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()
