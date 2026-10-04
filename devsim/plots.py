"""Generate report figures from results/*.csv into docs/figures/."""
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"figure.dpi": 150, "font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
C_N, C_P, C_OUT = "#1f6fd1", "#d1491f", "#2a9d5c"


def load(name):
    return np.genfromtxt(os.path.join(RES, name), delimiter=",", names=True)


def idvg():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for dev, col in (("nmos", C_N), ("pmos", C_P)):
        d = load(f"{dev}_idvg.csv")
        prm = json.load(open(os.path.join(RES, f"{dev}_params.json")))
        cols = d.dtype.names
        lab = "nMOS" if dev == "nmos" else "pMOS"
        for ax in axes:
            ax.plot(d[cols[0]], d[cols[2]] * 1e6, col, lw=2, label=f"{lab} |Vd| = 1.8 V")
            ax.plot(d[cols[0]], d[cols[1]] * 1e6, col, ls="--", lw=1.5, label=f"{lab} |Vd| = 0.05 V")
        axes[1].text(0.98, 0.32 - (0 if dev == "nmos" else 0.16),
                     f"{lab}: Vt={prm['Vt_maxgm']:.2f} V, SS={prm['SS_mV_per_dec']:.0f} mV/dec\n"
                     f"   DIBL={prm['DIBL_mV_per_V']:.0f} mV/V, Ion={prm['Ion_uA_per_um']:.0f} uA/um",
                     transform=axes[1].transAxes, va="top", ha="right", fontsize=8, color=col)
    axes[0].set(xlabel="|Vg| (V)", ylabel="|Id| (uA/um)", title="Transfer characteristics (linear)")
    axes[1].set(xlabel="|Vg| (V)", ylabel="|Id| (uA/um)", yscale="log", title="Transfer characteristics (log)",
                ylim=(1e-6, 2e3))
    axes[0].legend(fontsize=8)
    fig.suptitle("180 nm nMOS / pMOS  (L = 180 nm, tox = 4 nm)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "idvg.png"))


def idvd():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, dev, col in ((axes[0], "nmos", C_N), (axes[1], "pmos", C_P)):
        d = load(f"{dev}_idvd.csv")
        cols = d.dtype.names
        for k, c in enumerate(cols[1:]):
            vg = c.replace("Id_Vg", "")
            ax.plot(d[cols[0]], d[c] * 1e6, color=col, alpha=0.35 + 0.13 * k, lw=2, label=f"|Vg| = {vg} V")
        ax.set(xlabel="|Vd| (V)", title=("nMOS" if dev == "nmos" else "pMOS") + " output characteristics")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("|Id| (uA/um)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "idvd.png"))


def vtc():
    files = sorted(glob.glob(os.path.join(RES, "inverter_vtc_wp*.csv")))
    if not files:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    main = os.path.join(RES, "inverter_vtc_wp1.5.csv")
    for f in files:
        wp = os.path.basename(f)[len("inverter_vtc_wp"):-4]
        d = load(os.path.basename(f))
        is_main = f == main
        axes[0].plot(d["Vin"], d["Vout"], lw=2.5 if is_main else 1.2, alpha=1 if is_main else 0.6,
                     label=f"Wp/Wn = {wp}/1")
    axes[0].plot([0, 1.8], [0, 1.8], "k:", lw=1)
    d = load("inverter_vtc_wp1.5.csv")
    m = json.load(open(os.path.join(RES, "inverter_dc_wp1.5.json")))
    for x, lab in ((m["VIL"], "VIL"), (m["VM"], "VM"), (m["VIH"], "VIH")):
        axes[0].axvline(x, color="gray", ls="--", lw=0.8)
        axes[0].text(x, 1.86, lab, ha="center", fontsize=8)
    axes[0].set(xlabel="Vin (V)", ylabel="Vout (V)", xlim=(0, 1.8), ylim=(-0.05, 1.95),
                title="Voltage transfer characteristic")
    axes[0].legend(fontsize=8, loc="lower left")
    axes[0].text(1.0, 1.45, f"VM = {m['VM']:.3f} V\nNML = {m['NML']:.3f} V\nNMH = {m['NMH']:.3f} V\n"
                 f"gain = {m['peak_gain']:.1f}", fontsize=8, bbox=dict(fc="white", ec="0.8"))
    gain = -np.gradient(d["Vout"], d["Vin"])
    axes[1].plot(d["Vin"], gain, C_OUT, lw=2, label="|dVout/dVin|")
    axes[1].set(xlabel="Vin (V)", ylabel="|gain|", title="Gain and supply current (Wp/Wn = 1.5)")
    ax2 = axes[1].twinx()
    ax2.plot(d["Vin"], np.abs(d["Idd_A"]) * 1e6, color="#7a3cc4", lw=1.5, ls="--", label="IDD")
    ax2.set_ylabel("IDD (uA)")
    ax2.grid(False)
    h1, l1 = axes[1].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axes[1].legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    fig.suptitle("CMOS inverter, mixed-mode TCAD (VDD = 1.8 V)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "vtc.png"))


def tran():
    f = os.path.join(RES, "inverter_tran_wp1.5.csv")
    if not os.path.exists(f):
        return
    d = load("inverter_tran_wp1.5.csv")
    m = json.load(open(os.path.join(RES, "inverter_tran_wp1.5.json")))
    t = d["t_s"] * 1e12
    fig, ax = plt.subplots(figsize=(9, 3.8))
    ax.plot(t, d["Vin"], color="0.4", lw=1.5, label="Vin")
    ax.plot(t, d["Vout"], C_OUT, lw=2.2, label="Vout")
    ax.axhline(0.9, color="gray", ls=":", lw=0.8)
    ax.set(xlabel="time (ps)", ylabel="V", title=f"Transient response, CL = {m['CL_fF']:.0f} fF, "
           f"Wn/Wp = {m['Wn_um']:g}/{m['Wp_um']:g} um")
    ax.text(0.99, 0.5, f"tpHL = {m['tpHL_ps']:.1f} ps\ntpLH = {m['tpLH_ps']:.1f} ps\n"
            f"tp = {m['tp_avg_ps']:.1f} ps\ntf = {m['tfall_10_90_ps']:.1f} ps\ntr = {m['trise_10_90_ps']:.1f} ps",
            transform=ax.transAxes, ha="right", va="center", fontsize=9, bbox=dict(fc="white", ec="0.8"))
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "transient.png"))


if __name__ == "__main__":
    idvg()
    idvd()
    vtc()
    tran()
    print("figures in", FIG)
