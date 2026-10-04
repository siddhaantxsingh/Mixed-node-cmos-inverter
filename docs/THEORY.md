# Theory and design notes (viva guide)

## 1. What "mixed-mode" simulation means

There are three levels at which an inverter can be simulated:

| Level | Transistor described by | Accuracy | Speed |
|---|---|---|---|
| SPICE | Compact model equations (BSIM etc.) with fitted parameters | Only as good as the fit | ms |
| **Mixed-mode TCAD** | **Full 2D/3D device physics** (Poisson + continuity equations on a mesh) | Physical: geometry and doping directly set the behaviour | minutes–hours |
| Full TCAD of the whole circuit | Everything meshed, interconnect included | Highest | impractical |

In **mixed mode**, each transistor is a meshed TCAD device, and SDevice adds the **circuit equations** (Kirchhoff's current law at every node, plus the source and capacitor equations) to the same Newton system. At every bias point it solves *simultaneously*:

- **Poisson's equation** in each device: ∇·(ε∇ψ) = −q(p − n + N_D − N_A)
- **Electron and hole continuity**: ∇·J_n = q(R − G) + q ∂n/∂t, and the same for holes
- **KCL** at circuit nodes: Σ (device terminal currents + source currents + C dV/dt) = 0

The terminal current of a device (integrated over its contact) appears in the KCL of the node it connects to, and the node voltage becomes the boundary condition of that contact. That two-way coupling is what the `System { }` block in `inverter_*_des.cmd` sets up:

```
NMOS nmos1 ( "drain"=out "gate"=in "source"=0  "body"=0 )
PMOS pmos1 ( "drain"=out "gate"=in "source"=dd "body"=dd )
```

**Why use it?** You get inverter behaviour (VTC, delay) straight from the device design. Change Xj, tox, doping or workfunction and see the circuit effect with no compact model in between. It is how new technologies are evaluated before any SPICE model exists.

## 2. The CMOS inverter

- **Vin = 0**: nMOS off, pMOS on, so Vout = VDD (pulled up).
- **Vin = VDD**: nMOS on, pMOS off, so Vout = 0 (pulled down).
- In both states one transistor is off, so **static current is only leakage** (pA). That is the reason CMOS dominates.
- **During switching** both conduct briefly: the short-circuit current peak (118 µA in the reference run) occurs near VM.

### VTC regions (Vin rising)

| Region | nMOS | pMOS | Vout |
|---|---|---|---|
| A: Vin < Vtn | off | linear | VDD |
| B | saturation | linear | slightly below VDD |
| C: Vin ≈ VM | saturation | saturation | falls steeply (high gain) |
| D | linear | saturation | slightly above 0 |
| E: Vin > VDD − \|Vtp\| | linear | off | 0 |

### Key metrics (how `analyze_plt.py` computes them)

- **VM (switching threshold)**: the point where Vout = Vin.
- **VIL, VIH**: the two points where dVout/dVin = −1 (unity gain).
- **VOH = Vout(VIL)** and **VOL = Vout(VIH)**.
- **Noise margins**: NML = VIL − VOL and NMH = VOH − VIH. These are how much noise an input can carry while the next gate still reads it correctly.
- **Gain**: max |dVout/dVin| in region C. It is finite because of channel-length modulation and DIBL (a finite output resistance).
- **tpHL**: time from Vin crossing 50% (rising) to Vout crossing 50% (falling). **tpLH** is the reverse. tp = (tpHL + tpLH)/2.
- **Rise/fall time**: 10%–90% transition time of Vout.

### Switching threshold and sizing

With both devices in saturation at VM and Idn = Idp:

VM ≈ (Vtn + r·(VDD − |Vtp|)) / (1 + r),  where r = √(k_p / k_n) = √(μp·Wp / μn·Wn)

For VM = VDD/2 with |Vtn| = |Vtp|, you need r = 1, i.e. **matched drive currents**. Holes are slower than electrons, so the pMOS is made wider. Here the pMOS saturation current is 362 µA/µm against 529 µA/µm for the nMOS, so **Wp/Wn ≈ 529/362 ≈ 1.46**. Using 1.5 gives VM = 0.891 V, essentially VDD/2 = 0.9 V, and balanced noise margins (0.605 / 0.636 V).

The textbook ratio is about 2–3 (from μn/μp ≈ 2.5–3). In a 180 nm device it's smaller because **velocity saturation** limits electrons more than holes: in short channels, current is set more by v_sat than by low-field mobility. You can see this in the simulation.

### Why Vout overshoots at the start of each edge

In the transient plot, Vout rises slightly *above* VDD (peak 1.849 V, i.e. +49 mV) just as Vin starts rising, and dips to −37 mV on the falling input edge. When the input rises, the gate-drain overlap and fringe capacitance (C_gd) couples the input edge into the output node before the nMOS has turned on enough to pull it down. That is **Miller feedthrough**. Mixed-mode TCAD captures it automatically because the gate-drain coupling comes from the 2D electrostatics, not from a capacitance someone typed in.

### Delay estimate (sanity check)

tp ≈ CL·(VDD/2) / I_avg. With CL = 10 fF, VDD/2 = 0.9 V and an average discharge current of roughly 330 µA (between the nMOS saturation current and its linear-region current at Vout = 0.9 V), tp ≈ 27 ps. The simulation gives 27.5 ps.

## 3. Device design choices

| Choice | Reason |
|---|---|
| **L = 180 nm** | Task requirement; typical of the 0.18 µm node |
| **tox = 4 nm** | Typical 180 nm-node gate oxide (EOT ~3.5–4 nm). Thinner oxide gives more gate control, less DIBL and more drive |
| **VDD = 1.8 V** | Standard 180 nm supply |
| **Body 3×10¹⁷ cm⁻³** | Sets Vt and controls short-channel effects. The manual's 1×10¹⁷ is for a teaching device; at L = 180 nm it gives too much DIBL and leakage |
| **Workfunctions 4.25 / 4.97 eV** | Near band-edge metals, like n+/p+ poly (4.17 / 5.17 eV). Tuned so \|Vt\| ≈ 0.42 V for both devices, giving a symmetric inverter |
| **Gaussian S/D, Xj = 120 nm** | Same profile type as the manual. Deep, heavily doped S/D give low series resistance |
| **Extensions (LDD) 5×10¹⁸, Xj = 35 nm** | Shallow junctions next to the channel reduce charge sharing (short-channel effect) and the peak drain field (hot carriers) |
| **Si₃N₄ spacers, 100 nm** | Self-align the deep S/D away from the gate edge; the extensions sit under the spacers |
| **TiN gate (via workfunction)** | As in the manual; the metal is represented by the gate contact with `Workfunction=` |
| **AreaFactor** | 2D simulation gives current per µm of width; AreaFactor = W (µm) scales it to a real device |

### Workfunction to threshold voltage (sanity check)

For the nMOS: Φms = Φm − (χ + Eg/2 + φF) = 4.25 − (4.05 + 0.56 + 0.45) ≈ −0.81 V.

With φF = (kT/q)·ln(N_A/nᵢ) = 0.445 V, Cox = ε_ox/tox = 8.6×10⁻⁷ F/cm² and Qdep = √(2qε_si·N_A·2φF) = 3.0×10⁻⁷ C/cm²:

Vt = Φms + 2φF + Qdep/Cox = −0.81 + 0.89 + 0.35 ≈ **0.43 V**, against the simulated 0.42 V (max-gm). The long-channel formula works well at L = 180 nm because DIBL is small. The pMOS mirrors this with Φm = 4.97 eV.

## 4. Physics models (and why)

| Model (SDevice keyword) | Purpose |
|---|---|
| `EffectiveIntrinsicDensity(OldSlotboom)` | Bandgap narrowing in heavily doped S/D |
| `Mobility(DopingDep)` | Ionized-impurity scattering lowers mobility in doped regions |
| `HighFieldsaturation(GradQuasiFermi)` | Velocity saturation (~10⁷ cm/s). Essential at L = 180 nm; it is why Id grows roughly linearly with Vg instead of quadratically |
| `Enormal` | Surface-roughness scattering from the vertical gate field; reduces channel mobility at high Vg |
| `SRH(DopingDep TempDependence)` | Generation/recombination, sets junction leakage |

## 5. Reading the results

- **Subthreshold swing 75 mV/dec** (ideal limit 60): SS = 60·(1 + Cdep/Cox). 75 mV/dec means good gate control.
- **DIBL ~25 mV/V**: small. The drain barely lowers the source barrier at L = 180 nm with these extensions.
- **Ion/Ioff ~5×10⁷**: excellent standby behaviour.
- **Peak gain ~19**: set by output resistance (CLM + DIBL). A higher gain gives a sharper VTC and larger noise margins.
- **NML ≈ NMH**: confirms the sizing is balanced.

## 6. Exercises worth adding to the report

1. **Sizing sweep**: Wp/Wn = 1, 1.5, 2. VM moves up with Wp (`inverter.py dc --wp 2`).
2. **Supply scaling**: VDD = 1.2 V. Delay rises sharply as VDD approaches Vt.
3. **Load**: CL = 5, 10, 20 fF. tp scales almost linearly with CL (tp ≈ 0.69·R_eq·CL).
4. **Oxide thickness**: 3 vs 5 nm. Vt, SS, DIBL and delay change.
5. **Channel doping**: 1×10¹⁷ (manual value) vs 3×10¹⁷. More leakage and DIBL at low doping.

## 7. Likely viva questions

1. *What is mixed-mode simulation and how does it differ from SPICE?* (§1)
2. *Why is the pMOS wider?* Hole drive is lower; match drive currents for VM = VDD/2 (§2).
3. *Why is the ratio 1.5 and not 2.5?* Velocity saturation narrows the n/p current gap at short L.
4. *Define the noise margins. Where do VIL and VIH come from?* Unity-gain points (§2).
5. *Why is static power ~0 in CMOS?* One device is always off.
6. *Where does current peak in the VTC, and why?* At VM, where both devices are on (short-circuit current).
7. *What sets Vt here?* Gate workfunction, body doping, tox (§3 calculation).
8. *Why the LDD extensions and spacers?* Short-channel and hot-carrier control.
9. *What does AreaFactor do?* Converts per-µm 2D current to a W-µm device.
10. *Why does tpHL differ from tpLH?* Unequal pull-down and pull-up strength, even after sizing.
