*******************************************************************
*  NMOS (L = 180 nm) Id-Vd characteristics  --  SDevice
*  Run:  sdevice nmos_idvd_des.cmd      (needs nmos_180_msh.tdr from SDE)
*  Same structure as Lab Manual Exp 4/5 decks; VDD = 1.8 V.
*******************************************************************
File {
  Grid    = "nmos_180_msh.tdr"
  Plot    = "idvd_nmos_180"
  Current = "idvd_nmos_180"
  Output  = "idvd_nmos_180"
}

Electrode {
  { Name="gate"   Voltage= 0.0  Workfunction= 4.25 }
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
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson }
  Coupled ( Iterations= 100 LineSearchDamping= 1e-8 ){ Poisson electron hole }
  Save ( FilePrefix= "nmos_eq" )

  Load ( FilePrefix= "nmos_eq" )
  NewCurrentPrefix= "tmp_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Name= "gate" Voltage= 0.6 }
  ){ Coupled { Poisson electron hole } }
  NewCurrentPrefix= "IdVd_Vg0.6_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.05
    Goal { Name= "drain" Voltage= 1.8 }
  ){ Coupled { Poisson electron hole } }

  Load ( FilePrefix= "nmos_eq" )
  NewCurrentPrefix= "tmp_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Name= "gate" Voltage= 0.9 }
  ){ Coupled { Poisson electron hole } }
  NewCurrentPrefix= "IdVd_Vg0.9_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.05
    Goal { Name= "drain" Voltage= 1.8 }
  ){ Coupled { Poisson electron hole } }

  Load ( FilePrefix= "nmos_eq" )
  NewCurrentPrefix= "tmp_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Name= "gate" Voltage= 1.2 }
  ){ Coupled { Poisson electron hole } }
  NewCurrentPrefix= "IdVd_Vg1.2_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.05
    Goal { Name= "drain" Voltage= 1.8 }
  ){ Coupled { Poisson electron hole } }

  Load ( FilePrefix= "nmos_eq" )
  NewCurrentPrefix= "tmp_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Name= "gate" Voltage= 1.5 }
  ){ Coupled { Poisson electron hole } }
  NewCurrentPrefix= "IdVd_Vg1.5_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.05
    Goal { Name= "drain" Voltage= 1.8 }
  ){ Coupled { Poisson electron hole } }

  Load ( FilePrefix= "nmos_eq" )
  NewCurrentPrefix= "tmp_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.1
    Goal { Name= "gate" Voltage= 1.8 }
  ){ Coupled { Poisson electron hole } }
  NewCurrentPrefix= "IdVd_Vg1.8_"
  Quasistationary(
    InitialStep= 1e-3 Increment= 1.5 MinStep= 1e-6 MaxStep= 0.05
    Goal { Name= "drain" Voltage= 1.8 }
  ){ Coupled { Poisson electron hole } }
}
