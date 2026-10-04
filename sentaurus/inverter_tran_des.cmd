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
*  Run:  sdevice inverter_tran_des.cmd
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
  Output = "inverter_tran"
}

System {
  Vsource_pset vdd ( dd 0 ) { dc = 0.0 }
  *- input pulse: 0 -> 1.8 V, delay 200 ps, rise/fall 20 ps, high for 400 ps
  *  pulse = ( v_low  v_high  delay  t_rise  t_fall  t_high  period )
  Vsource_pset vin ( in 0 ) { pulse = ( 0.0 1.8 200e-12 20e-12 20e-12 400e-12 2e-9 ) }
  Capacitor_pset cl ( out 0 ) { capacitance = 10e-15 }      * 10 fF load
  NMOS nmos1 ( "drain"=out "gate"=in "source"=0  "body"=0 )
  PMOS pmos1 ( "drain"=out "gate"=in "source"=dd "body"=dd )
  Plot "inverter_tran_sys.plt" ( time() v(in) v(out) v(dd) i(vdd dd) i(cl out) )
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
  Transient= BE
}

Solve {
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson }
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson Electron Hole Contact Circuit }

  *- DC operating point at VDD = 1.8 V, Vin = 0 (Vout = 1.8 V)
  Quasistationary(
    InitialStep= 1e-2 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Parameter= vdd.dc Voltage= 1.8 }
  ){ Coupled { Poisson Electron Hole Contact Circuit } }

  *- transient: 0 .. 1 ns, step <= 2 ps
  NewCurrentPrefix= "tran_"
  Transient(
    InitialTime= 0 FinalTime= 1.0e-9
    InitialStep= 1e-13 MaxStep= 2e-12 MinStep= 1e-17 Increment= 1.3
  ){ Coupled { Poisson Electron Hole Contact Circuit } }
}
