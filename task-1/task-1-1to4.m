(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 1                                                 *)
(* Short Calculations in Task 1.1, 1.2, and 1.4           *)
(* ====================================================== *)

(* Start FeynCalc *)
<< FeynCalc`

(* Task 1.1: Z-Boson Decay Trace *)
PolarizationVector[p, \[Mu]] . SpinorUBar[k1, 0] . 
 GA[\[Mu]] . (g - f*GA5) . SpinorV[k2, 0]
% ComplexConjugate[%];
DoPolarizationSums[%, p, ExtraFactor -> 1/3];
FermionSpinSum[%];
Simplify[DiracSimplify[%]]


(* Task 1.2: Photon Matrix Element *)
SpinorVBar[p2, 0] . GA[\[Mu]] . SpinorU[p1, 0] . SpinorUBar[p3, 0] . 
 GA[\[Mu]] . SpinorV[p4, 0]
% ComplexConjugate[%];
FermionSpinSum[%, ExtraFactor -> 1/4];
Simplify[DiracSimplify[%]]


(* Task 1.2: Z-Boson Matric Element *)
SpinorVBar[p2, 0] . GA[\[Mu]] . (a - b*GA5) . SpinorU[p1, 0] . 
 SpinorUBar[p3, 0] . GA[\[Mu]] . (a - b GA5) . SpinorV[p4, 0]
% ComplexConjugate[%];
FermionSpinSum[%];
Simplify[DiracSimplify[%]]


(* Task 1.4: t-Wb Trace *)
FCClearScalarProducts[]
SP[p2] = \!\(\*SubsuperscriptBox[\(m\), \(W\), \(2\)]\);

SpinorUBar[p3, Subscript[m, b]] . GA[\[Mu]] . (1 - GA5) . 
 SpinorU[p1, Subscript[m, t]] . 
 Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p2, -I]]]
% ComplexConjugate[%];
DoPolarizationSums[%, p2];
FermionSpinSum[%];
Simplify[DiracSimplify[%]]


(* Task 1.4: t-Hb Trace *)
SpinorUBar[p3, m3] . (1 + GA5) . SpinorU[p1, m1]
% ComplexConjugate[%];
FermionSpinSum[%, ExtraFactor -> 1/2];
Simplify[DiracSimplify[%]]
