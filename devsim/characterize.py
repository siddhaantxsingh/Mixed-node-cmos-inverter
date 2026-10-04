"""
Single-device characterisation of the 180 nm nMOS and pMOS
(equivalent to the lab's Exp 4/5 IdVg / IdVd SDevice runs).

Outputs (results/):
  <dev>_idvg.csv     Vg, Id at |Vd| = 0.05 V and 1.8 V
  <dev>_idvd.csv     Vd, Id for |Vg| = 0.6 ... 1.8 V
  <dev>_params.json  Vt (max-gm & constant-current), SS, DIBL, Ion, Ioff
"""

import json
import os
import sys
import time

import numpy as np

os.environ.setdefault("DEVSIM_MATH_LIBS", "/usr/local/lib/libmkl_rt.so.3")
import devsim as ds  # noqa: E402

import tcad  # noqa: E402

VDD = 1.8
OUT = os.path.join(os.path.dirname(__file__), "..", "results")
A_PER_UM = 1e-4  # DEVSIM 2D current is A/cm of width


def fresh(name, params):
    for d in ds.get_device_list():
        ds.delete_device(device=d)
    for m in ds.get_mesh_list():
        ds.delete_mesh(mesh=m)
    tcad.create_device(name, params)
    tcad.solve()
    tcad.enable_dd(name, params)
    tcad.solve()


def idvg(name, params, vd, s):
    """Sweep |Vg| 0 -> VDD at fixed |Vd|; returns arrays (Vg, |Id|)."""
    fresh(name, params)
    tcad.ramp_param(name, f"{name}_drain_bias", s * vd, step=0.1)
    vgs, ids = [0.0], [abs(tcad.contact_current(name, "drain")) * A_PER_UM]
    target = s * VDD
    vg = 0.0
    while abs(vg) < VDD - 1e-9:
        nxt = vg + s * 0.025
        tcad.ramp_param(name, f"{name}_gate_bias", nxt, step=0.025)
        vg = nxt
        vgs.append(abs(vg))
        ids.append(abs(tcad.contact_current(name, "drain")) * A_PER_UM)
    return np.array(vgs), np.array(ids)


def idvd(name, params, vg_list, s):
    curves = {}
    for vg in vg_list:
        fresh(name, params)
        tcad.ramp_param(name, f"{name}_gate_bias", s * vg, step=0.1)
        vds, ids = [0.0], [0.0]
        vd = 0.0
        while abs(vd) < VDD - 1e-9:
            nxt = vd + s * 0.05
            tcad.ramp_param(name, f"{name}_drain_bias", nxt, step=0.05)
            vd = nxt
            vds.append(abs(vd))
            ids.append(abs(tcad.contact_current(name, "drain")) * A_PER_UM)
        curves[vg] = (np.array(vds), np.array(ids))
    return curves


def extract(vg, id_lin, id_sat, vd_lin=0.05):
    """Threshold voltage, subthreshold swing, DIBL, Ion, Ioff."""
    gm = np.gradient(id_lin, vg)
    k = int(np.argmax(gm))
    vt_gm = vg[k] - id_lin[k] / gm[k] - vd_lin / 2  # max-gm linear extrapolation
    # constant-current Vt at 100 nA/um (W/L-normalised, common for 180 nm)
    icc = 1e-7 * 1.0

    def vt_cc(i):
        return float(np.interp(np.log10(icc), np.log10(np.maximum(i, 1e-30)), vg))

    vt_lin, vt_sat = vt_cc(id_lin), vt_cc(id_sat)
    # subthreshold swing: min of dVg/dlog10(Id) in the subthreshold region
    lg = np.log10(np.maximum(id_lin, 1e-30))
    mask = (id_lin > 1e-12) & (id_lin < icc)
    ss = np.gradient(vg, lg)[mask]
    ss_mv = float(np.min(ss[ss > 0]) * 1000) if np.any(ss > 0) else float("nan")
    return dict(
        Vt_maxgm=float(vt_gm),
        Vt_cc_lin=vt_lin,
        Vt_cc_sat=vt_sat,
        DIBL_mV_per_V=float((vt_lin - vt_sat) / (VDD - vd_lin) * 1000),
        SS_mV_per_dec=ss_mv,
        Ion_uA_per_um=float(id_sat[-1] * 1e6),
        Ioff_pA_per_um=float(id_sat[0] * 1e12),
        Ion_Ioff=float(id_sat[-1] / max(id_sat[0], 1e-30)),
    )


def run(kind, **overrides):
    name = "nmos" if kind == "n" else "pmos"
    params = tcad.MOSParams(kind, **overrides)
    s = 1.0 if kind == "n" else -1.0
    t0 = time.time()
    vg, il = idvg(name, params, 0.05, s)
    _, isat = idvg(name, params, VDD, s)
    np.savetxt(os.path.join(OUT, f"{name}_idvg.csv"), np.c_[vg, il, isat],
               delimiter=",", header="Vg_abs,Id_lin_A_per_um(Vd=0.05),Id_sat_A_per_um(Vd=1.8)", comments="")
    curves = idvd(name, params, [0.6, 0.9, 1.2, 1.5, 1.8], s)
    vd = curves[1.8][0]
    np.savetxt(os.path.join(OUT, f"{name}_idvd.csv"), np.c_[[vd] + [c[1] for c in curves.values()]].T,
               delimiter=",", header="Vd_abs," + ",".join(f"Id_Vg{v}" for v in curves), comments="")
    res = extract(vg, il, isat)
    res["workfunction_eV"] = params.wf
    res["body_doping_cm3"] = params.body
    res["tox_nm"] = params.tox * 1000
    res["L_nm"] = params.L * 1000
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, f"{name}_params.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(name, json.dumps(res, indent=2))
    return res


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    kinds = sys.argv[1:] or ["n", "p"]
    for k in kinds:
        run(k)
