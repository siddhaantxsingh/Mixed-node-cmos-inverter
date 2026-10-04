*******************************************************************
*  CMOS INVERTER -- MIXED-MODE (device + circuit) simulation, SDevice
*  nMOS (L=180 nm, W=1.0 um) + pMOS (L=180 nm, W=1.5 um), VDD = 1.8 V
*
*        vdd ----+---------------
*                |  pMOS  (source,body = vdd)
*        in  ----+-- gates        out --[ CL ]-- 0      (CL: transient only)
*                |  nMOS  (source,body = 0)
*        0   ----+---------------
*
*  Both transistors are full 2D TCAD devices (meshes from SDE); the
*  System{} block wires them into a circuit and SDevice solves the device
*  equations and the circuit (KCL) equations together.
*
*  Run:  sdevice inverter_dc_des.cmd
*******************************************************************
Device NMOS {
  Electrode {
    { Name="source" Voltage= 0.0 }
    { Name="drain"  Voltage= 0.0 }
    { Name="gate"   Voltage= 0.0  Workfunction= 4.25 }
    { Name="body"   Voltage= 0.0 }
  }
  File {
    Grid = "nmos_180_msh.tdr"
    Plot = "inv_nmos"
  }
  Physics {
    AreaFactor= 1.0     * device width in um (2D current is per um)
    EffectiveIntrinsicDensity( OldSlotboom )
    Mobility(
      DopingDep
      eHighFieldsaturation( GradQuasiFermi )
      hHighFieldsaturation( GradQuasiFermi )
      Enormal
    )
    Recombination( SRH( DopingDep TempDependence ) )
  }
  Plot {
    *--Density and Currents, etc
    eDensity hDensity
    TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
    eMobility hMobility
    eVelocity hVelocity
    eQuasiFermi hQuasiFermi
    *--Fields and charges
    ElectricField/Vector Potential SpaceCharge
    *--Doping Profiles
    Doping DonorConcentration AcceptorConcentration
    *--Generation/Recombination
    SRH Band2BandGeneration * Auger
    *--Driving forces
    eGradQuasiFermi/Vector hGradQuasiFermi/Vector
    eEparallel hEparallel eENormal hENormal
    *--Band structure/Composition
    BandGap BandGapNarrowing Affinity
    ConductionBand ValenceBand
  }
}

Device PMOS {
  Electrode {
    { Name="source" Voltage= 0.0 }
    { Name="drain"  Voltage= 0.0 }
    { Name="gate"   Voltage= 0.0  Workfunction= 4.97 }
    { Name="body"   Voltage= 0.0 }
  }
  File {
    Grid = "pmos_180_msh.tdr"
    Plot = "inv_pmos"
  }
  Physics {
    AreaFactor= 1.5     * device width in um (2D current is per um)
    EffectiveIntrinsicDensity( OldSlotboom )
    Mobility(
      DopingDep
      eHighFieldsaturation( GradQuasiFermi )
      hHighFieldsaturation( GradQuasiFermi )
      Enormal
    )
    Recombination( SRH( DopingDep TempDependence ) )
  }
  Plot {
    *--Density and Currents, etc
    eDensity hDensity
    TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
    eMobility hMobility
    eVelocity hVelocity
    eQuasiFermi hQuasiFermi
    *--Fields and charges
    ElectricField/Vector Potential SpaceCharge
    *--Doping Profiles
    Doping DonorConcentration AcceptorConcentration
    *--Generation/Recombination
    SRH Band2BandGeneration * Auger
    *--Driving forces
    eGradQuasiFermi/Vector hGradQuasiFermi/Vector
    eEparallel hEparallel eENormal hENormal
    *--Band structure/Composition
    BandGap BandGapNarrowing Affinity
    ConductionBand ValenceBand
  }
}

File {
  Output = "inverter_dc"
}

System {
  Vsource_pset vdd ( dd 0 ) { dc = 0.0 }
  Vsource_pset vin ( in 0 ) { dc = 0.0 }
  NMOS nmos1 ( "drain"=out "gate"=in "source"=0  "body"=0 )
  PMOS pmos1 ( "drain"=out "gate"=in "source"=dd "body"=dd )
  Plot "inverter_dc_sys.plt" ( v(in) v(out) v(dd) i(vdd dd) i(nmos1 out) i(pmos1 out) )
}

Math {
  Extrapolate
  Derivatives
  RelErrControl
  Digits=5
  ErrRef(electron)=1.e10
  ErrRef(hole)=1.e10
  Iterations=20
  Notdamped=100
  Method=Blocked
  SubMethod=Super
  ACMethod=Blocked
  ACSubMethod=Super
}

Solve {
  *- equilibrium (all sources at 0 V)
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson }
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson Electron Hole Contact Circuit }

  *- ramp the supply to VDD = 1.8 V (Vin = 0, so Vout follows VDD)
  Quasistationary(
    InitialStep= 1e-2 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Parameter= vdd.dc Voltage= 1.8 }
  ){ Coupled { Poisson Electron Hole Contact Circuit } }

  *- VTC: sweep Vin 0 -> 1.8 V in <= 10 mV steps
  NewCurrentPrefix= "vtc_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.3 MinStep= 1e-7 MaxStep= 0.01
    Goal { Parameter= vin.dc Voltage= 1.8 }
  ){ Coupled { Poisson Electron Hole Contact Circuit }
     Plot ( FilePrefix= "vtc_vin" Time= ( 0.5 ) NoOverwrite )
  }
}
