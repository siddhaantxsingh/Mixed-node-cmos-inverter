"""
2D drift-diffusion MOSFET simulator on DEVSIM, mirroring the Sentaurus
SDE/SDevice setup of the SDPS lab (Exp 4/5) with L = 180 nm.

Physics (same model set as the lab's SDevice decks, simplified):
  * Poisson + electron/hole continuity, Scharfetter-Gummel discretisation
  * Doping-dependent mobility (Caughey-Thomas / Masetti-like)
  * High-field velocity saturation (Canali form, edge-parallel field)
  * SRH recombination
  * Metal gate (TiN) with workfunction, ideal ohmic S/D/body contacts
  * Optional circuit coupling for mixed-mode (contacts attach to circuit nodes)

Units: DEVSIM works in cm, V, A, s. Device currents are per cm of width
(2D); scripts convert to A/um.
"""

import math
import os

os.environ.setdefault("DEVSIM_MATH_LIBS", "/usr/local/lib/libmkl_rt.so.3")
import devsim as ds  # noqa: E402

# Use UMFPACK instead of MKL PARDISO: with circuit nodes in the system the
# matrix mixes O(1) circuit entries with O(1e-14) device entries, and
# PARDISO's pivoting reported spurious zero pivots (see docs/devsim-notes.md).
import devsim.umfpack.umfshim  # noqa: E402,F401

# ----------------------------------------------------------------- constants
q = 1.602176634e-19
k_B = 1.380649e-23
eps0 = 8.8541878128e-14  # F/cm
T = 300.0
Vt = k_B * T / q
ni = 1.0e10  # cm^-3 (intrinsic carrier density used throughout)
EPS = {"Silicon": 11.7, "Oxide": 3.9, "Nitride": 7.5}
CHI_SI = 4.05  # electron affinity (eV)
EG_SI = 1.12
UM = 1e-4  # um -> cm


class MOSParams:
    """Geometry/doping of one transistor (lengths in um, doping in cm^-3)."""

    def __init__(self, kind, L=0.18, tox=0.004, spacer=0.10, sd_len=0.35,
                 depth=2.0, body=3e17, sd_peak=5e19, sd_junction_depth=0.12,
                 ext_peak=5e18, ext_junction_depth=0.035, lateral_factor=0.8,
                 workfunction=None, width_factor=1.0):
        assert kind in ("n", "p")
        self.kind = kind
        self.L, self.tox, self.spacer, self.sd_len, self.depth = L, tox, spacer, sd_len, depth
        self.body = body
        self.sd_peak, self.sd_jd = sd_peak, sd_junction_depth
        self.ext_peak, self.ext_jd = ext_peak, ext_junction_depth
        self.lat = lateral_factor
        # TiN-like gate workfunctions tuned per device type (see docs).
        self.wf = workfunction if workfunction is not None else (4.25 if kind == "n" else 4.97)
        self.width_factor = width_factor  # W relative to 1 um (scales currents & charges)

    # x coordinates (um) of feature edges, left to right
    @property
    def xs(self):
        x0 = 0.0
        x1 = x0 + self.sd_len          # source / spacer edge
        x2 = x1 + self.spacer          # gate left edge
        x3 = x2 + self.L               # gate right edge
        x4 = x3 + self.spacer          # spacer / drain edge
        x5 = x4 + self.sd_len
        return x0, x1, x2, x3, x4, x5


# NOTE: DEVSIM's parser binds unary minus tighter than '^', so
# exp(-(a)^2) means exp((-a)^2). Always write exp(0 - (a)^2).
def _sigma(peak, junction, depth):
    """Gaussian std-dev such that peak*exp(-(d/s)^2) = junction at d = depth."""
    return depth / math.sqrt(math.log(peak / junction))


def build_mesh(name, p: MOSParams):
    """Tensor-product mesh: silicon body, gate oxide, nitride spacers."""
    x0, x1, x2, x3, x4, x5 = p.xs
    ds.create_2d_mesh(mesh=name)

    # --- x lines: fine near junctions & channel
    def xl(pos, ps, ns=None):
        ds.add_2d_mesh_line(mesh=name, dir="x", pos=pos * UM, ps=ps * UM, ns=(ns or ps) * UM)

    xl(x0, 0.05)
    xl(x1 - 0.05, 0.02)
    xl(x1, 0.006)
    xl(x2, 0.003)
    xl((x2 + x3) / 2, 0.01)
    xl(x3, 0.003)
    xl(x4, 0.006)
    xl(x4 + 0.05, 0.02)
    xl(x5, 0.05)

    # --- y lines: y = 0 is the Si surface, +y into silicon, -y above
    def yl(pos, ps, ns=None):
        ds.add_2d_mesh_line(mesh=name, dir="y", pos=pos * UM, ps=ps * UM, ns=(ns or ps) * UM)

    sp_top = -max(p.tox, 0.010)
    yl(sp_top, 0.002)
    yl(-p.tox, 0.001)
    yl(0.0, 0.0005)
    yl(0.01, 0.002)
    yl(0.05, 0.005)
    yl(0.2, 0.02)
    yl(0.5, 0.08)
    yl(p.depth, 0.3)
    yl(p.depth + 0.01, 0.01)

    um = lambda v: v * UM
    # "gas" regions surround the device; DEVSIM contacts lie on the boundary
    # between a simulated region and gas. The gas over the gate oxide stands
    # in for the TiN gate electrode (not simulated, represented by a contact).
    ds.add_2d_region(mesh=name, material="gas", region="gas_top", xl=um(x0), xh=um(x5), yl=um(sp_top), yh=0.0)
    ds.add_2d_region(mesh=name, material="gas", region="gas_bot", xl=um(x0), xh=um(x5), yl=um(p.depth), yh=um(p.depth + 0.01))
    ds.add_2d_region(mesh=name, material="Silicon", region="si",
                     xl=um(x0), xh=um(x5), yl=0.0, yh=um(p.depth))
    ds.add_2d_region(mesh=name, material="Oxide", region="ox",
                     xl=um(x2), xh=um(x3), yl=um(-p.tox), yh=0.0)
    ds.add_2d_region(mesh=name, material="Nitride", region="spL",
                     xl=um(x1), xh=um(x2), yl=um(sp_top), yh=0.0)
    ds.add_2d_region(mesh=name, material="Nitride", region="spR",
                     xl=um(x3), xh=um(x4), yl=um(sp_top), yh=0.0)

    ds.add_2d_interface(mesh=name, name="si_ox", region0="si", region1="ox",
                        xl=um(x2), xh=um(x3), yl=0.0, yh=0.0, bloat=1e-10)
    ds.add_2d_interface(mesh=name, name="si_spL", region0="si", region1="spL",
                        xl=um(x1), xh=um(x2), yl=0.0, yh=0.0, bloat=1e-10)
    ds.add_2d_interface(mesh=name, name="si_spR", region0="si", region1="spR",
                        xl=um(x3), xh=um(x4), yl=0.0, yh=0.0, bloat=1e-10)
    # No oxide/spacer side-wall interfaces: the triple points (Si/ox/nitride)
    # would get redundant continuity equations and a singular Jacobian.
    # The side walls are Neumann boundaries instead (negligible effect).

    # S/D contacts stop 2 nm short of the spacers so no node is both a
    # contact node and a Si/nitride interface node.
    ds.add_2d_contact(mesh=name, name="source", material="metal", region="si",
                      xl=um(x0), xh=um(x1 - 0.002), yl=0.0, yh=0.0, bloat=1e-10)
    ds.add_2d_contact(mesh=name, name="drain", material="metal", region="si",
                      xl=um(x4 + 0.002), xh=um(x5), yl=0.0, yh=0.0, bloat=1e-10)
    ds.add_2d_contact(mesh=name, name="body", material="metal", region="si",
                      xl=um(x0), xh=um(x5), yl=um(p.depth), yh=um(p.depth), bloat=1e-10)
    ds.add_2d_contact(mesh=name, name="gate", material="metal", region="ox",
                      xl=um(x2), xh=um(x3), yl=um(-p.tox), yh=um(-p.tox), bloat=1e-10)
    ds.finalize_mesh(mesh=name)
    ds.create_device(mesh=name, device=name)


# ----------------------------------------------------------------- helpers
def regions(dev):
    """Simulated regions (excludes the surrounding gas)."""
    return [r for r in ds.get_region_list(device=dev) if ds.get_material(device=dev, region=r) != "gas"]


def node_model(dev, reg, name, expr, *deriv_vars):
    ds.node_model(device=dev, region=reg, name=name, equation=expr)
    for v in deriv_vars:
        ds.node_model(device=dev, region=reg, name=f"{name}:{v}", equation=f"diff({expr},{v})")


def edge_model(dev, reg, name, expr, *deriv_vars):
    ds.edge_model(device=dev, region=reg, name=name, equation=expr)
    for v in deriv_vars:
        for n in ("n0", "n1"):
            ds.edge_model(device=dev, region=reg, name=f"{name}:{v}@{n}", equation=f"diff({expr},{v}@{n})")


def edge_from_node(dev, reg, model):
    ds.edge_from_node_model(device=dev, region=reg, node_model=model)


def set_scaling(dev, region, k):
    """Scale the 'third dimension' (device width) by k: multiplies every flux
    and every stored charge in the region, i.e. W -> k*W."""
    ds.edge_model(device=dev, region=region, name="ScaledEdgeCouple", equation=f"{k}*EdgeCouple")
    ds.node_model(device=dev, region=region, name="ScaledNodeVolume", equation=f"{k}*NodeVolume")
    ds.set_parameter(device=dev, region=region, name="edge_couple_model", value="ScaledEdgeCouple")
    ds.set_parameter(device=dev, region=region, name="node_volume_model", value="ScaledNodeVolume")


# ----------------------------------------------------------------- doping
def doping(dev, p: MOSParams):
    x0, x1, x2, x3, x4, x5 = [v * UM for v in p.xs]
    s_sd = _sigma(p.sd_peak, p.body, p.sd_jd) * UM
    s_ext = _sigma(p.ext_peak, p.body, p.ext_jd) * UM

    def window(xl, xh, s_lat):
        # Gaussian lateral decay outside [xl, xh] (Sentaurus "Gaussian Factor")
        return (f"ifelse(x < {xl}, exp(0 - ((x-({xl}))/{s_lat})^2),"
                f" ifelse(x > {xh}, exp(0 - ((x-({xh}))/{s_lat})^2), 1))")

    def gauss(peak, s, xl, xh):
        return f"{peak}*exp(0 - (y/{s})^2)*{window(xl, xh, p.lat * s)}"

    sd = f"({gauss(p.sd_peak, s_sd, x0, x1)} + {gauss(p.sd_peak, s_sd, x4, x5)})"
    ext = f"({gauss(p.ext_peak, s_ext, x1, x2)} + {gauss(p.ext_peak, s_ext, x3, x4)})"
    implant = f"({sd} + {ext})"
    if p.kind == "n":
        donors, acceptors = implant, f"{p.body}"
    else:
        donors, acceptors = f"{p.body}", implant
    node_model(dev, "si", "Donors", donors)
    node_model(dev, "si", "Acceptors", acceptors)
    node_model(dev, "si", "NetDoping", "Donors - Acceptors")
    node_model(dev, "si", "TotalDoping", "Donors + Acceptors")


# ----------------------------------------------------------------- physics
# DEVSIM sign conventions (match devsim.python_packages.simple_physics):
#   Poisson:     div(eps*E) + q*(n - p - N) = 0      -> node model +q(n-p-N)
#   Electrons:   div(Jn) - q*U = q dn/dt             -> generation -qU, charge -qn
#   Holes:       div(Jp) + q*U = -q dp/dt            -> generation +qU, charge +qp
#   SG:  Jn = q*mu*Vt/L*(n1*B(-v) - n0*B(v)),  Jp = -q*mu*Vt/L*(p1*B(v) - p0*B(-v)),
#        v = (psi0 - psi1)/Vt
def set_params(dev):
    for reg in regions(dev):
        mat = ds.get_material(device=dev, region=reg)
        ds.set_parameter(device=dev, region=reg, name="Permittivity", value=EPS[mat] * eps0)
    ds.set_parameter(device=dev, name="q", value=q)
    ds.set_parameter(device=dev, name="V_t", value=Vt)
    ds.set_parameter(device=dev, name="n_i", value=ni)
    # SRH lifetimes
    ds.set_parameter(device=dev, name="taun", value=1e-7)
    ds.set_parameter(device=dev, name="taup", value=1e-7)
    # Caughey-Thomas doping dependence (Si, 300 K, Masetti-like constants)
    for k, v in dict(mun_min=52.2, mun_max=1417.0, Nrefn=9.68e16, alphan=0.68,
                     mup_min=44.9, mup_max=470.5, Nrefp=2.23e17, alphap=0.719,
                     vsatn=1.07e7, vsatp=8.37e6, betan=1.109, betap=1.213).items():
        ds.set_parameter(device=dev, name=k, value=v)


def potential_only(dev):
    """Equilibrium Poisson in all regions (used for the initial guess)."""
    for reg in regions(dev):
        ds.node_solution(device=dev, region=reg, name="Potential")
        ds.edge_from_node_model(device=dev, region=reg, node_model="Potential")
        edge_model(dev, reg, "ElectricField", "(Potential@n0 - Potential@n1)*EdgeInverseLength", "Potential")
        edge_model(dev, reg, "DField", "Permittivity*ElectricField", "Potential")
    # Initial guess: charge-neutral potential in silicon (Newton diverges
    # from psi = 0 with 5e19 doping); oxide/spacers start at 0.
    node_model(dev, "si", "NeutralPotential", "V_t*asinh(NetDoping/(2*n_i))")
    ds.set_node_values(device=dev, region="si", name="Potential", init_from="NeutralPotential")
    # Silicon: Boltzmann carriers from potential (equilibrium)
    node_model(dev, "si", "IntrinsicElectrons", "n_i*exp(Potential/V_t)", "Potential")
    node_model(dev, "si", "IntrinsicHoles", "n_i^2/IntrinsicElectrons", "Potential")
    node_model(dev, "si", "PotentialIntrinsicCharge",
               "q*(IntrinsicElectrons - IntrinsicHoles - NetDoping)", "Potential")
    for reg in regions(dev):
        if reg == "si":
            ds.equation(device=dev, region=reg, name="PotentialEquation", variable_name="Potential",
                        node_model="PotentialIntrinsicCharge", edge_model="DField",
                        variable_update="log_damp")
        else:
            ds.equation(device=dev, region=reg, name="PotentialEquation", variable_name="Potential",
                        edge_model="DField", variable_update="log_damp")
    for itf in ds.get_interface_list(device=dev):
        name = f"continuous_potential_{itf}"
        ds.interface_model(device=dev, interface=itf, name=name, equation="Potential@r0-Potential@r1")
        ds.interface_model(device=dev, interface=itf, name=f"{name}:Potential@r0", equation="1")
        ds.interface_model(device=dev, interface=itf, name=f"{name}:Potential@r1", equation="-1")
        ds.interface_equation(device=dev, interface=itf, name="PotentialEquation",
                              interface_model=name, type="continuous")


def bias_name(dev, contact):
    return f"{dev}_{contact}_bias"


def contacts(dev, p: MOSParams, circuit_nodes=None, drift_diffusion=False):
    """Ohmic S/D/body and metal gate contacts.

    circuit_nodes: optional dict contact -> circuit node name (mixed mode).
    Without it, each contact's voltage is the device parameter
    `<dev>_<contact>_bias`.
    """
    circuit_nodes = circuit_nodes or {}
    # Equilibrium carrier densities at ohmic contacts, branch-free:
    #   n = ni*exp(asinh(N/2ni)) = (N + sqrt(N^2+4ni^2))/2 ,  p = ni*exp(-asinh(N/2ni))
    # exact, no cancellation for either sign of N, and no ifelse (DEVSIM
    # evaluated an ifelse model lazily inside contact assembly with the wrong
    # branch; see docs/devsim-notes.md).
    node_model(dev, "si", "asinhN", "asinh(NetDoping/(2*n_i))")
    node_model(dev, "si", "celec", "n_i*exp(asinhN)")
    node_model(dev, "si", "chole", "n_i*exp(0 - asinhN)")
    for m in ("asinhN", "celec", "chole"):  # force evaluation now
        ds.get_node_model_values(device=dev, region="si", name=m)
    phi_ms_offset = p.wf - (CHI_SI + EG_SI / 2)  # gate workfunction vs intrinsic level

    for c in ("source", "drain", "body", "gate"):
        node = circuit_nodes.get(c)
        bias = node if node else bias_name(dev, c)
        if not node:
            ds.set_parameter(device=dev, name=bias, value=0.0)
        region = "ox" if c == "gate" else "si"
        if c == "gate":
            expr = f"Potential - ({bias}) + {phi_ms_offset}"
        else:
            expr = f"Potential - ({bias}) - V_t*asinhN"
        cm = f"{c}_pot"
        ds.contact_node_model(device=dev, contact=c, name=cm, equation=expr)
        ds.contact_node_model(device=dev, contact=c, name=f"{cm}:Potential", equation="1")
        kw = {}
        if node:
            ds.contact_node_model(device=dev, contact=c, name=f"{cm}:{node}", equation="-1")
            kw["circuit_node"] = node
        ds.contact_equation(device=dev, contact=c, name="PotentialEquation", node_model=cm,
                            edge_charge_model="DField", **kw)
        if drift_diffusion and c != "gate":
            for var, eq, cur, eqm in (
                ("Electrons", "ElectronContinuityEquation", "ElectronCurrent",
                 "Electrons - celec"),
                ("Holes", "HoleContinuityEquation", "HoleCurrent",
                 "Holes - chole"),
            ):
                nm = f"{c}_{var}"
                ds.contact_node_model(device=dev, contact=c, name=nm, equation=eqm)
                ds.contact_node_model(device=dev, contact=c, name=f"{nm}:{var}", equation="1")
                ds.contact_equation(device=dev, contact=c, name=eq, node_model=nm,
                                    edge_current_model=cur, **kw)


def drift_diffusion(dev):
    """Add electron/hole continuity in silicon (after an equilibrium solve)."""
    reg = "si"
    for v in ("Electrons", "Holes"):
        ds.node_solution(device=dev, region=reg, name=v)
        ds.edge_from_node_model(device=dev, region=reg, node_model=v)
    ds.set_node_values(device=dev, region=reg, name="Electrons", init_from="IntrinsicElectrons")
    ds.set_node_values(device=dev, region=reg, name="Holes", init_from="IntrinsicHoles")

    # Poisson with free carriers
    node_model(dev, reg, "PotentialNodeCharge", "q*(Electrons - Holes - NetDoping)", "Electrons", "Holes")
    ds.equation(device=dev, region=reg, name="PotentialEquation", variable_name="Potential",
                node_model="PotentialNodeCharge", edge_model="DField", variable_update="log_damp")

    # SRH
    srh = "(Electrons*Holes - n_i^2)/(taup*(Electrons + n_i) + taun*(Holes + n_i))"
    node_model(dev, reg, "USRH", srh, "Electrons", "Holes")
    node_model(dev, reg, "ElectronGeneration", "-q*USRH", "Electrons", "Holes")
    node_model(dev, reg, "HoleGeneration", "q*USRH", "Electrons", "Holes")
    node_model(dev, reg, "NCharge", "-q*Electrons", "Electrons")
    node_model(dev, reg, "PCharge", "q*Holes", "Holes")

    # Doping-dependent low-field mobility (node), averaged onto edges
    node_model(dev, reg, "mu_n_dop", "mun_min + (mun_max - mun_min)/(1 + (TotalDoping/Nrefn)^alphan)")
    node_model(dev, reg, "mu_p_dop", "mup_min + (mup_max - mup_min)/(1 + (TotalDoping/Nrefp)^alphap)")
    edge_from_node(dev, reg, "mu_n_dop")
    edge_from_node(dev, reg, "mu_p_dop")
    ds.edge_model(device=dev, region=reg, name="mu_n0", equation="0.5*(mu_n_dop@n0 + mu_n_dop@n1)")
    ds.edge_model(device=dev, region=reg, name="mu_p0", equation="0.5*(mu_p_dop@n0 + mu_p_dop@n1)")

    # Velocity saturation (Canali), driven by edge field magnitude.
    # sqrt(E^2 + 1) avoids a non-differentiable abs() at E = 0.
    efield = "((Potential@n0 - Potential@n1)*EdgeInverseLength)"
    emag = f"(({efield})^2 + 1.0)^0.5"
    mun = f"(mu_n0/(1 + (mu_n0*{emag}/vsatn)^betan)^(1/betan))"
    mup = f"(mu_p0/(1 + (mu_p0*{emag}/vsatp)^betap)^(1/betap))"

    # Scharfetter-Gummel currents
    # Bernoulli terms as named edge models with explicit derivatives:
    # DEVSIM's symbolic diff() does not differentiate B() inline.
    ds.edge_model(device=dev, region=reg, name="vdiff", equation="(Potential@n0 - Potential@n1)/V_t")
    ds.edge_model(device=dev, region=reg, name="vdiff:Potential@n0", equation="V_t^(-1)")
    ds.edge_model(device=dev, region=reg, name="vdiff:Potential@n1", equation="-V_t^(-1)")
    ds.edge_model(device=dev, region=reg, name="Bern01", equation="B(vdiff)")
    ds.edge_model(device=dev, region=reg, name="Bern01:Potential@n0", equation="dBdx(vdiff)*vdiff:Potential@n0")
    ds.edge_model(device=dev, region=reg, name="Bern01:Potential@n1", equation="-Bern01:Potential@n0")
    ds.edge_model(device=dev, region=reg, name="Bern10", equation="Bern01 + vdiff")      # = B(-vdiff)
    ds.edge_model(device=dev, region=reg, name="Bern10:Potential@n0", equation="Bern01:Potential@n0 + vdiff:Potential@n0")
    ds.edge_model(device=dev, region=reg, name="Bern10:Potential@n1", equation="Bern01:Potential@n1 + vdiff:Potential@n1")
    # kahan3 = compensated 3-term sum: avoids cancellation in n1*B(-v) - n0*B(v)
    # when n ~ 5e19 (otherwise Newton stalls at ~1e-5 relative error).
    jn = f"q*{mun}*EdgeInverseLength*V_t*kahan3(Electrons@n1*Bern01, Electrons@n1*vdiff, -Electrons@n0*Bern01)"
    jp = f"-q*{mup}*EdgeInverseLength*V_t*kahan3(Holes@n1*Bern01, -Holes@n0*Bern01, -Holes@n0*vdiff)"
    edge_model(dev, reg, "ElectronCurrent", jn, "Potential", "Electrons", "Holes")
    edge_model(dev, reg, "HoleCurrent", jp, "Potential", "Electrons", "Holes")

    ds.equation(device=dev, region=reg, name="ElectronContinuityEquation", variable_name="Electrons",
                time_node_model="NCharge", edge_model="ElectronCurrent",
                node_model="ElectronGeneration", variable_update="positive")
    ds.equation(device=dev, region=reg, name="HoleContinuityEquation", variable_name="Holes",
                time_node_model="PCharge", edge_model="HoleCurrent",
                node_model="HoleGeneration", variable_update="positive")


def contact_current(dev, contact):
    """Total terminal current (A per cm width, scaled by width_factor)."""
    if contact == "gate":
        return 0.0
    return (ds.get_contact_current(device=dev, contact=contact, equation="ElectronContinuityEquation")
            + ds.get_contact_current(device=dev, contact=contact, equation="HoleContinuityEquation"))


def create_device(name, p: MOSParams, circuit_nodes=None):
    """Build mesh, doping, physics; solve equilibrium; enable drift-diffusion."""
    build_mesh(name, p)
    if p.width_factor != 1.0:
        for reg in regions(name):
            set_scaling(name, reg, p.width_factor)
    set_params(name)
    doping(name, p)
    potential_only(name)
    contacts(name, p, circuit_nodes=circuit_nodes, drift_diffusion=False)
    return name


def enable_dd(name, p, circuit_nodes=None):
    drift_diffusion(name)
    contacts(name, p, circuit_nodes=circuit_nodes, drift_diffusion=True)


def solve(abs_err=1e10, rel_err=1e-9, max_iter=40, **kw):
    ds.solve(type="dc", absolute_error=abs_err, relative_error=rel_err, maximum_iterations=max_iter, **kw)


def ramp_param(dev, name, target, step=0.05, min_step=1e-4, callback=None):
    """Ramp a bias parameter with adaptive step (halve on failure)."""
    v = ds.get_parameter(device=dev, name=name)
    while abs(target - v) > 1e-12:
        dv = math.copysign(min(step, abs(target - v)), target - v)
        try:
            ds.set_parameter(device=dev, name=name, value=v + dv)
            solve()
            v += dv
            if callback:
                callback(v)
            step = min(step * 1.5, 0.1)
        except ds.error:
            ds.set_parameter(device=dev, name=name, value=v)
            step /= 2
            if step < min_step:
                raise RuntimeError(f"ramp of {name} failed at {v}")
    return v
