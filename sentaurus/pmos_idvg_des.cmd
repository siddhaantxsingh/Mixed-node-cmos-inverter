*******************************************************************
*  PMOS (L = 180 nm) Id-Vg characteristics  --  SDevice
*  Run:  sdevice pmos_idvg_des.cmd      (needs pmos_180_msh.tdr from SDE)
*  Same structure as Lab Manual Exp 4/5 decks; VDD = 1.8 V.
*******************************************************************
File {
  Grid    = "pmos_180_msh.tdr"
  Plot    = "idvg_pmos_180"
  Current = "idvg_pmos_180"
  Output  = "idvg_pmos_180"
}

Electrode {
  { Name="gate"   Voltage= 0.0  Workfunction= 4.97 }
  { Name="drain"  Voltage= 0.0 }
  { Name="source" Voltage= 0.0 }
  { Name="body"   Voltage= 0.0 }
}

Physics {
  EffectiveIntrinsicDensity( OldSlotboom )
  Mobility(
    DopingDep
    eHighFieldsaturation( GradQuasiFermi )
    hHighFieldsaturation( GradQuasiFermi )
    Enormal
  )
  Recombination(
    SRH( DopingDep TempDependence )
  )
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
  *- initial solution
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson }
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson electron hole }

  *- Vd = -0.05 V (linear region)
  Quasistationary(
    DoZero
    InitialStep= 1e-3 Increment= 1.5
    MinStep= 1e-6 MaxStep= 0.025
    Goal { Name= "drain" Voltage= -0.05 }
  ){ Coupled { Poisson electron hole } }
  Save ( FilePrefix= "pmos_vd_lin" )

  NewCurrentPrefix= "IdVg_lin_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5
    MinStep= 1e-6 MaxStep= 0.025
    Goal { Name= "gate" Voltage= -1.8 }
  ){ Coupled { Poisson electron hole } }

  *- Vd = VDD (saturation)
  Load ( FilePrefix= "pmos_vd_lin" )
  NewCurrentPrefix= "tmp_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5
    MinStep= 1e-6 MaxStep= 0.1
    Goal { Name= "drain" Voltage= -1.8 }
  ){ Coupled { Poisson electron hole } }
  NewCurrentPrefix= "IdVg_sat_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5
    MinStep= 1e-6 MaxStep= 0.025
    Goal { Name= "gate" Voltage= -1.8 }
  ){ Coupled { Poisson electron hole } }
}
