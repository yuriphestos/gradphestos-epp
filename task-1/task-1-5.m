(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 1                                                 *)
(* Task 1.5: Decay Width of 125 GeV SM Higgs              *)
(* ====================================================== *)

(* Start FeynCalc *)
<< FeynCalc`

(* Partial width of h-boson to ZZ* *)

Clear[mZ]
Clear[mh]
Clear[Msq]


(* final state Z polarization matrix, summed over polarizations *)
Pmunu := 3*1/
   3 (FourVector[p3, \[Nu]]*FourVector[p3, \[Mu]]/mZ^2 - 
    MetricTensor[\[Mu], \[Nu]])
Pmunu


(*Z and fermions vertex*)
ZF[\[Mu]_] := (1/2) GA[\[Mu]] . (gv - ga GA5)


(* Final state Fermionic tensor *)
FT := Tr[GS[p1] . ZF[\[Mu]] . GS[p2] . ZF[\[Nu]]]
FT


(*"1"- fermion, "2" - fermion, "3" - Z. \
Three body matrix element squared, without coupling constants, \
that will be taken care of later. Quantity in the denominator - \
Z propagator *)

1/(2 Pair[Momentum[p1], Momentum[p2]] - mZ^2)^2 Contract[FT*Pmunu]



(* Matrix elements squared (coupling constants are dealt with separate\
ly). Kinematic substitutions are done to make \
|M|^2 to depend only on m_{12}^2 and m_{23}^2 *)
Msq = 1/(2 Pair[Momentum[p1], Momentum[p2]] - mZ^2)^2 Contract[
     FT*Pmunu] /. {Pair[Momentum[p2], Momentum[p2]] -> 0,
    Pair[Momentum[p1], Momentum[p1]] -> 0, 
    Pair[Momentum[p3], Momentum[p3]] -> mZ^2,
    Pair[Momentum[p1], Momentum[p2]] -> m12sq/2, 
    Pair[Momentum[p1], 
      Momentum[p3]] -> (mh^2 + mZ^2 - m23sq - m12sq)/2 - mZ^2/2, 
    Pair[Momentum[p2], Momentum[p3]] -> m23sq/2 - mZ^2/2} // 
  FullSimplify



(* Integration over m23sq and m12sq *)
$Assumptions = 
 mh > 0 && mZ > 0 && m12sq < mh^2 && mh^2 - m12sq - mZ^2 > 0 && 
  m12sq > 0



(* Limits of integration over m_{23}^2 at fixed m_{12}^2 *)

E2st := (m12sq)/(2*Sqrt[m12sq])
E3st := (mh^2 - m12sq - mZ^2)/(2*Sqrt[m12sq])
m23sqMax := (E2st + E3st)^2 - (Sqrt[E2st^2] - Sqrt[E3st^2 - mZ^2])^2
m23sqMin := (E2st + E3st)^2 - (Sqrt[E2st^2] + Sqrt[E3st^2 - mZ^2])^2
m23sqMin // FullSimplify
m23sqMax  // FullSimplify


(* Matrix element integrated over m_{23}^2 *)

Msq1 = Integrate[Msq, {m23sq, m23sqMin, m23sqMax}] // FullSimplify


(* Matrix element integrated over m_{23}^2 and m_{12}^2 *)
mH := 125.
mZ := 91.18
MsqNumeric = 
 NIntegrate[
  Msq1 /. {mh -> 125, mZ -> 91.18, gv -> -0.038, ga -> 0.5}, {m12sq, 
   0., (mH - mZ)^2}]


(* Number of fermionic decay channels, one pair from on-shell Z, \
the other from off-shell Z *)
Nchan = 2


(* Coupling constants. First multiplier - hZZ coupling squared, \
Second multiplier - Zff coupling *)
Couplings := (2 mZ^2/v)^2*(gW/cosw)^2 /. {gW^2 -> 4 Pi \[Alpha]W, 
   cosw -> 0.88}
Couplings


CouplingsN := Couplings /. {\[Alpha]W -> 1/137/0.223, v -> 246}
CouplingsN


(* Partial Decay Width, in GeV, using PDG formula 49.22 *)
GaZ = MsqNumeric * (1/(2*Pi)^3)/(32*(mH)^3)*Nchan*CouplingsN


(* Branching ratio of the decay Z to electons or muons *)
BrZll = 0.084 2/2.49


(* Higgs width in GeV *) 
GaH := 4.07*10^(-3)

(* Branching ratio of h -> 4l *)
Brh4l = GaZ*BrZll/GaH


(* SM Higgs branching ratio from CERN report 2016 *)
(* Link: \
https://twiki.cern.ch/twiki/bin/view/LHCPhysics/\
CERNYellowReportPageBR#H_llll_ll *)
Brh4lCERN := 1.24*10^(-4)
Brh4l/Brh4lCERN

GaZ*BrZll
