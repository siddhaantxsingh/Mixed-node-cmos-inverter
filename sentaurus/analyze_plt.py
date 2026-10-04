#!/usr/bin/env python3
"""
Extract CMOS-inverter metrics from Sentaurus SDevice mixed-mode output.

Usage (on the lab machine, after running the decks):
    python3 analyze_plt.py dc   inverter_dc_sys.plt
    python3 analyze_plt.py tran inverter_tran_sys.plt
    python3 analyze_plt.py list inverter_dc_sys.plt      # show dataset names

Reads the DF-ISE text .plt format written by SDevice (datasets [...] and
Data { ... }), so it needs only Python 3 (numpy optional, not required).

DC   -> VM, VIL, VIH, VOH, VOL, NML, NMH, peak gain, peak supply current
TRAN -> tpHL, tpLH, average tp, 10-90 % rise / fall times
"""
import re
import sys

VDD = 1.8


def read_plt(path):
    text = open(path, encoding="latin-1").read()
    m = re.search(r"datasets\s*=\s*\[(.*?)\]", text, re.S)
    if not m:
        sys.exit(f"{path}: no 'datasets' header found (is this a .plt file?)")
    names = re.findall(r'"([^"]*)"', m.group(1))
    d = re.search(r"Data\s*\{(.*)\}", text, re.S)
    nums = [float(x) for x in d.group(1).split()]
    n = len(names)
    rows = [nums[i:i + n] for i in range(0, len(nums) - n + 1, n)]
    return names, {name: [r[i] for r in rows] for i, name in enumerate(names)}


def pick(names, data, *patterns):
    """First dataset whose name matches all regex patterns (case-insensitive)."""
    for nm in names:
        if all(re.search(p, nm, re.I) for p in patterns):
            return data[nm], nm
    sys.exit(f"no dataset matching {patterns}; run 'list' to see names:\n  " + "\n  ".join(names))


def interp(x0, xs, ys):
    for i in range(1, len(xs)):
        a, b = xs[i - 1], xs[i]
        if (a - x0) * (b - x0) <= 0 and a != b:
            return ys[i - 1] + (x0 - a) / (b - a) * (ys[i] - ys[i - 1])
    return float("nan")


def dc(path):
    names, data = read_plt(path)
    vin, n_in = pick(names, data, r"^in\b|v\(in\)|\bin\s")
    vout, n_out = pick(names, data, r"^out\b|v\(out\)|\bout\s")
    print(f"using Vin = '{n_in}', Vout = '{n_out}'")
    # gain by central differences
    g = []
    for i in range(len(vin)):
        j0, j1 = max(i - 1, 0), min(i + 1, len(vin) - 1)
        g.append((vout[j1] - vout[j0]) / (vin[j1] - vin[j0]) if vin[j1] != vin[j0] else 0.0)
    diff = [o - i for i, o in zip(vin, vout)]
    vm = interp(0.0, diff, vin)
    idx = [i for i, x in enumerate(g) if x <= -1]
    vil = interp(-1.0, g[: idx[0] + 1], vin[: idx[0] + 1]) if idx else float("nan")
    vih = interp(-1.0, g[idx[-1]:], vin[idx[-1]:]) if idx else float("nan")
    voh, vol = interp(vil, vin, vout), interp(vih, vin, vout)
    res = dict(VM=vm, VIL=vil, VIH=vih, VOH=voh, VOL=vol, NML=vil - vol, NMH=voh - vih, peak_gain=-min(g))
    try:
        idd, n_i = pick(names, data, r"vdd")
        res["peak_Idd_uA"] = max(abs(x) for x in idd) * 1e6
    except SystemExit:
        pass
    for k, v in res.items():
        print(f"  {k:12s} {v:10.4f}")


def crossing(t, v, level, rising, after):
    for i in range(1, len(t)):
        if t[i] <= after:
            continue
        a, b = v[i - 1], v[i]
        if (rising and a < level <= b) or (not rising and a > level >= b):
            return t[i - 1] + (level - a) / (b - a) * (t[i] - t[i - 1])
    return float("nan")


def tran(path):
    names, data = read_plt(path)
    t, _ = pick(names, data, r"time")
    vin, _ = pick(names, data, r"^in\b|v\(in\)|\bin\s")
    vout, _ = pick(names, data, r"^out\b|v\(out\)|\bout\s")
    h = VDD / 2
    tin_r = crossing(t, vin, h, True, 0)
    tout_f = crossing(t, vout, h, False, tin_r - 1e-12)
    tin_f = crossing(t, vin, h, False, tin_r)
    tout_r = crossing(t, vout, h, True, tin_f - 1e-12)
    tf = crossing(t, vout, 0.1 * VDD, False, tin_r - 1e-12) - crossing(t, vout, 0.9 * VDD, False, tin_r - 1e-12)
    tr = crossing(t, vout, 0.9 * VDD, True, tin_f - 1e-12) - crossing(t, vout, 0.1 * VDD, True, tin_f - 1e-12)
    res = dict(tpHL_ps=(tout_f - tin_r) * 1e12, tpLH_ps=(tout_r - tin_f) * 1e12,
               tfall_ps=tf * 1e12, trise_ps=tr * 1e12)
    res["tp_avg_ps"] = (res["tpHL_ps"] + res["tpLH_ps"]) / 2
    for k, v in res.items():
        print(f"  {k:12s} {v:10.2f}")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("dc", "tran", "list"):
        sys.exit(__doc__)
    mode, path = sys.argv[1], sys.argv[2]
    if mode == "list":
        print("\n".join(read_plt(path)[0]))
    elif mode == "dc":
        dc(path)
    else:
        tran(path)
