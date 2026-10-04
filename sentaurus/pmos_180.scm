;; =====================================================================
;;  pMOS, L = 180 nm  --  Sentaurus Structure Editor (SDE) script
;;  SDPS Lab (ECE 3144) mini-project: mixed-mode CMOS inverter
;;
;;  Run:   sde -e -l pmos_180.scm
;;  Output: pmos_180_msh.tdr  (used by the SDevice decks)
;;
;;  Same construction as Lab Manual Exp 4/5 (nMOSFET / pMOSFET), scaled
;;  to a 180 nm channel. Coordinates in um; the silicon surface is at
;;  y = Dsub (as in the manual, where it was y = 5).
;; =====================================================================
(sde:clear)
(sdegeo:set-default-boolean "ABA")

;; ---------------- parameters ----------------------------------------
(define Lg    0.18)    ; gate (channel) length            180 nm
(define Tox   0.004)   ; gate oxide thickness             4 nm
(define Lsp   0.10)    ; Si3N4 spacer length              100 nm
(define Hsp   0.010)   ; spacer / gate-stack height       10 nm
(define Lsd   0.35)    ; source / drain contact length
(define Dsub  2.0)     ; silicon depth (body contact at y = 0)

(define Nbody   3e17)  ; PhosphorusActiveConcentration body doping
(define Nsd     5e19)  ; S/D peak (BoronActiveConcentration)
(define Xj_sd   0.12)  ; S/D junction depth
(define Next    5e18)  ; extension (LDD) peak
(define Xj_ext  0.035) ; extension junction depth
(define LatFac  0.8)   ; lateral Gaussian factor (manual: 0.8)

(define X0 0.0)
(define X1 (+ X0 Lsd))   ; source | spacer
(define X2 (+ X1 Lsp))   ; spacer | gate
(define X3 (+ X2 Lg))    ; gate   | spacer
(define X4 (+ X3 Lsp))   ; spacer | drain
(define X5 (+ X4 Lsd))
(define Ys Dsub)         ; silicon surface
(define Yox (+ Ys Tox))
(define Ytop (+ Ys Hsp))

;; ---------------- regions (names as in the manual) ------------------
(sdegeo:create-rectangle (position X0 0.0 0.0) (position X5 Ys 0.0)   "Silicon" "region_body")
(sdegeo:create-rectangle (position X2 Ys 0.0)  (position X3 Yox 0.0)  "SiO2"    "region_oxide")
(sdegeo:create-rectangle (position X2 Yox 0.0) (position X3 Ytop 0.0) "TiN"     "region_metal")
(sdegeo:create-rectangle (position X1 Ys 0.0)  (position X2 Ytop 0.0) "Si3N4"   "region_source")
(sdegeo:create-rectangle (position X3 Ys 0.0)  (position X4 Ytop 0.0) "Si3N4"   "region_drain")

;; ---------------- contacts ------------------------------------------
(sdegeo:define-contact-set "source" 4.0 (color:rgb 1.0 0.0 0.0) "##")
(sdegeo:define-contact-set "drain"  4.0 (color:rgb 0.0 1.0 0.0) "##")
(sdegeo:define-contact-set "gate"   4.0 (color:rgb 0.0 0.0 1.0) "##")
(sdegeo:define-contact-set "body"   4.0 (color:rgb 1.0 1.0 0.0) "##")

(sdegeo:set-current-contact-set "source")
(sdegeo:set-contact-edges (list (car (find-edge-id (position (/ (+ X0 X1) 2) Ys 0.0)))) "source")
(sdegeo:set-current-contact-set "drain")
(sdegeo:set-contact-edges (list (car (find-edge-id (position (/ (+ X4 X5) 2) Ys 0.0)))) "drain")
(sdegeo:set-current-contact-set "body")
(sdegeo:set-contact-edges (list (car (find-edge-id (position (/ X5 2) 0.0 0.0)))) "body")
;; Gate: the TiN region becomes the gate contact (its boundary), then the
;; metal body is removed -- SDevice models it through Workfunction.
(sdegeo:set-current-contact-set "gate")
(sdegeo:set-contact-boundary-edges (list (car (find-body-id (position (/ (+ X2 X3) 2) (/ (+ Yox Ytop) 2) 0.0)))) "gate")
(sdegeo:delete-region (list (car (find-body-id (position (/ (+ X2 X3) 2) (/ (+ Yox Ytop) 2) 0.0)))))

;; ---------------- doping --------------------------------------------
;; (a) constant body doping
(sdedr:define-constant-profile "ConstantProfileDefinition_body" "PhosphorusActiveConcentration" Nbody)
(sdedr:define-constant-profile-region "ConstantProfilePlacement_body" "ConstantProfileDefinition_body" "region_body")

;; (b) reference lines on the silicon surface (Mesh > Define Ref/Eval window > Line)
(sdedr:define-refeval-window "RefEvalWin_source"     "Line" (position X0 Ys 0.0) (position X1 Ys 0.0))
(sdedr:define-refeval-window "RefEvalWin_drain"      "Line" (position X4 Ys 0.0) (position X5 Ys 0.0))
(sdedr:define-refeval-window "RefEvalWin_source_ext" "Line" (position X1 Ys 0.0) (position X2 Ys 0.0))
(sdedr:define-refeval-window "RefEvalWin_drain_ext"  "Line" (position X3 Ys 0.0) (position X4 Ys 0.0))

;; (c) Gaussian profiles: peak at the surface, value Nbody at the junction depth
(sdedr:define-gaussian-profile "AnalyticalProfileDefinition_sd" "BoronActiveConcentration"
  "PeakPos" 0 "PeakVal" Nsd "ValueAtDepth" Nbody "Depth" Xj_sd "Gauss" "Factor" LatFac)
(sdedr:define-gaussian-profile "AnalyticalProfileDefinition_ext" "BoronActiveConcentration"
  "PeakPos" 0 "PeakVal" Next "ValueAtDepth" Nbody "Depth" Xj_ext "Gauss" "Factor" LatFac)

(sdedr:define-analytical-profile-placement "AnalyticalProfilePlacement_source"
  "AnalyticalProfileDefinition_sd" "RefEvalWin_source" "Both" "NoReplace" "Eval")
(sdedr:define-analytical-profile-placement "AnalyticalProfilePlacement_drain"
  "AnalyticalProfileDefinition_sd" "RefEvalWin_drain" "Both" "NoReplace" "Eval")
(sdedr:define-analytical-profile-placement "AnalyticalProfilePlacement_source_ext"
  "AnalyticalProfileDefinition_ext" "RefEvalWin_source_ext" "Both" "NoReplace" "Eval")
(sdedr:define-analytical-profile-placement "AnalyticalProfilePlacement_drain_ext"
  "AnalyticalProfileDefinition_ext" "RefEvalWin_drain_ext" "Both" "NoReplace" "Eval")

;; ---------------- meshing -------------------------------------------
;; (i) global
(sdedr:define-refeval-window "RefEvalWin_global" "Rectangle" (position X0 0.0 0.0) (position X5 Ytop 0.0))
(sdedr:define-refinement-size "RefinementDefinition_global" 0.05 0.2 0.01 0.01)
(sdedr:define-refinement-placement "RefinementPlacement_global" "RefinementDefinition_global" "RefEvalWin_global")
(sdedr:define-refinement-function "RefinementDefinition_global" "DopingConcentration" "MaxTransDiff" 1)
;; (ii) channel / junctions: fine mesh near the surface
(sdedr:define-refeval-window "RefEvalWin_channel" "Rectangle" (position (- X1 0.05) (- Ys 0.15) 0.0) (position (+ X4 0.05) Ys 0.0))
(sdedr:define-refinement-size "RefinementDefinition_channel" 0.005 0.005 0.001 0.0005)
(sdedr:define-refinement-placement "RefinementPlacement_channel" "RefinementDefinition_channel" "RefEvalWin_channel")
(sdedr:define-refinement-function "RefinementDefinition_channel" "MaxLenInt" "region_body" "region_oxide" 0.0005 1.4 "DoubleSide" "UseRegionNames")
(sdedr:define-refinement-function "RefinementDefinition_channel" "MaxLenInt" "region_body" "region_source" 0.001 1.4 "DoubleSide" "UseRegionNames")
(sdedr:define-refinement-function "RefinementDefinition_channel" "MaxLenInt" "region_body" "region_drain" 0.001 1.4 "DoubleSide" "UseRegionNames")

;; ---------------- save + build mesh ---------------------------------
(sde:save-model "pmos_180")
(sde:build-mesh "snmesh" "-a -c boxmethod" "pmos_180_msh")
