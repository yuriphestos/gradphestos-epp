(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 1                                                 *)
(* Task 1.6: e+e- -> W+W- Pair Production                  *)
(* ====================================================== *)


(* Start FeynCalc *)
<< FeynCalc`

(* e^+e^--> W^+W^- *)

(* Z-ff Vertex *)
ZF[\[Mu]_] := GA[\[Mu]] . (gv - ga*GA5)

ZF[\[Mu]]

(* Triple Gauge Vertex *)
TGV[\[Mu]_, \[Nu]_, \[Lambda]_, p3_, p4_, q1_] := 
 MetricTensor[\[Mu], \[Nu]] . FV[p3 - p4, \[Lambda]] + 
  MetricTensor[\[Lambda], \[Nu]] . FV[-q1 - p3, \[Mu]] + 
  MetricTensor[\[Lambda], \[Mu]] . FV[q1 + p4, \[Nu]]

TGV[\[Mu], \[Nu], \[Lambda], p3, p4, q1]


q1 = p1 + p2;
q2 = p1 - p3;


(* Couplings *)
M\[Gamma]cp = -e^2/ScalarProduct[q1, q1]
MZcp = gW/(2*cosw)*gWWZ/(ScalarProduct[q1, q1] - mZ^2)
M\[Nu]cp = -gW^2/(4*ScalarProduct[q2, q2]) 


(* M\[Gamma]^2 *)
FCClearScalarProducts[]
SP[p3, p3] = mW^2;
SP[p4, p4] = mW^2;

SpinorVBar[p2, 0] . GA[\[Lambda]] . SpinorU[p1, 0] . 
  TGV[\[Mu], \[Nu], \[Lambda], p3, p4, q1] . 
  Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p4, -I]]] . 
  Pair[LorentzIndex[\[Nu]], Momentum[Polarization[p3, -I]]];
% ComplexConjugate[%];
DoPolarizationSums[%, p3];
DoPolarizationSums[%, p4];
FermionSpinSum[%, ExtraFactor -> 1/4];

M\[Gamma]sq = M\[Gamma]cp^2*Simplify[DiracSimplify[%]]


(* MZ^2 *)
SpinorVBar[p2, 0] . ZF[\[Lambda]] . SpinorU[p1, 0] . 
  TGV[\[Mu], \[Nu], \[Lambda], p3, p4, q1] . 
  Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p4, -I]]] . 
  Pair[LorentzIndex[\[Nu]], Momentum[Polarization[p3, -I]]];
% ComplexConjugate[%];
DoPolarizationSums[%, p3];
DoPolarizationSums[%, p4];
FermionSpinSum[%, ExtraFactor -> 1/4];

MZsq = MZcp^2*Simplify[DiracSimplify[%]]


(* M\[Nu]^2 *)
SpinorVBar[p2, 0] . GA[\[Mu]] . GS[q2] . GA[\[Nu]] . (1 - GA5) . 
  SpinorU[p1, 0] . Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p4, -I]]] . 
  Pair[LorentzIndex[\[Nu]], Momentum[Polarization[p3, -I]]];
% ComplexConjugate[%];
FermionSpinSum[%, ExtraFactor -> 1/4];
DoPolarizationSums[%, p3];
DoPolarizationSums[%, p4];

M\[Nu]sq = M\[Nu]cp^2*Simplify[DiracSimplify[%]]


(* Subscript[M, \[Gamma]]Subscript[M, Z] *)
ComplexConjugate[
   SpinorVBar[p2, 0] . GA[\[Lambda]] . SpinorU[p1, 0] . 
    TGV[\[Mu], \[Nu], \[Lambda], p3, p4, q1] . 
    Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p4, -I]]] . 
    Pair[LorentzIndex[\[Nu]], Momentum[Polarization[p3, -I]]]] . 
  SpinorVBar[p2, 0] . ZF[\[Sigma]] . SpinorU[p1, 0] . 
  TGV[\[Alpha], \[Beta], \[Sigma], p3, p4, q1] . 
  Pair[LorentzIndex[\[Alpha]], Momentum[Polarization[p4, -I]]] . 
  Pair[LorentzIndex[\[Beta]], Momentum[Polarization[p3, -I]]];
DoPolarizationSums[%, p3];
DoPolarizationSums[%, p4];
FermionSpinSum[%, ExtraFactor -> 1/4];

M\[Gamma]MZ = M\[Gamma]cp*MZcp*Simplify[DiracSimplify[%]]


(* Subscript[M, Z]Subscript[M, \[Nu]] *)
ComplexConjugate[
   SpinorVBar[p2, 0] . ZF[\[Lambda]] . SpinorU[p1, 0] . 
    TGV[\[Mu], \[Nu], \[Lambda], p3, p4, q1] . 
    Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p4, -I]]] . 
    Pair[LorentzIndex[\[Nu]], Momentum[Polarization[p3, -I]]]] . 
  SpinorVBar[p2, 0] . GA[\[Alpha]] . GS[q2] . GA[\[Beta]] . (1 - GA5) . 
  SpinorU[p1, 0] . 
  Pair[LorentzIndex[\[Alpha]], Momentum[Polarization[p4, -I]]] . 
  Pair[LorentzIndex[\[Beta]], Momentum[Polarization[p3, -I]]];
DoPolarizationSums[%, p3];
DoPolarizationSums[%, p4];
FermionSpinSum[%, ExtraFactor -> 1/4];

MZM\[Nu] = M\[Nu]cp*MZcp*Simplify[DiracSimplify[%]]


(* Subscript[M, \[Gamma]]Subscript[M, \[Nu]] *)
ComplexConjugate[
   SpinorVBar[p2, 0] . GA[\[Lambda]] . SpinorU[p1, 0] . 
    TGV[\[Mu], \[Nu], \[Lambda], p3, p4, q1] . 
    Pair[LorentzIndex[\[Mu]], Momentum[Polarization[p4, -I]]] . 
    Pair[LorentzIndex[\[Nu]], Momentum[Polarization[p3, -I]]]] . 
  SpinorVBar[p2, 0] . GA[\[Alpha]] . GS[q2] . GA[\[Beta]] . (1 - GA5) . 
  SpinorU[p1, 0] . 
  Pair[LorentzIndex[\[Alpha]], Momentum[Polarization[p4, -I]]] . 
  Pair[LorentzIndex[\[Beta]], Momentum[Polarization[p3, -I]]];
DoPolarizationSums[%, p3];
DoPolarizationSums[%, p4];
FermionSpinSum[%, ExtraFactor -> 1/4];

M\[Gamma]M\[Nu] = M\[Gamma]cp *M\[Nu]cp*Simplify[DiracSimplify[%]]


(* Defining Scalar Product of Momenta *)
p = Sqrt[Ecm^2/4 - mW^2];
ScalarProduct[p1, p1] = 0;
ScalarProduct[p2, p2] = 0;
ScalarProduct[p3, p3] = mW^2;
ScalarProduct[p4, p4] = mW^2;

ScalarProduct[p1, p2] = Ecm^2/2;
ScalarProduct[p1, p3] = Ecm/2*(Ecm/2 - p*Cos[\[Theta]]);
ScalarProduct[p2, p4] = Ecm/2*(Ecm/2 - p*Cos[\[Theta]]);
ScalarProduct[p1, p4] = Ecm/2*(Ecm/2 + p*Cos[\[Theta]]);
ScalarProduct[p2, p3] = Ecm/2*(Ecm/2 + p*Cos[\[Theta]]);
ScalarProduct[p3, p4] = Ecm^2/2 - mW^2;

Pair[Momentum[p1 + p2], Momentum[p1 + p2]] = Ecm^2;
Pair[Momentum[p1 - p3], Momentum[p1 - p3]] = 
  mW^2 - 2*Ecm/2*(Ecm/2 - p*Cos[\[Theta]]);
Eps[Momentum[p1], Momentum[p2], Momentum[p3], Momentum[p4]] = 0;

(* Total Amplitude Squared (Only \[Gamma],Z exchanges considered) *)
Msq = M\[Gamma]sq + MZsq + 2*(M\[Gamma]MZ) // FullSimplify


(* Msq with numerical values *)
mZ = 91.18;
mW = 80.37;
gv = -0.038;
ga = -0.5;
gW = 0.63;
e = Sqrt[4*Pi/137];
cosw = 0.88;

Msq


(* Total Cross-section in unit pb *)
MInt = Integrate[Msq*Sin[\[Theta]], {\[Theta], 0, \[Pi]}] // FullSimplify

Subscript[\[Sigma], total] = (p*(MInt))/(16*Pi*Ecm^3)*3.89379*10^8


(* Final plot with gWWZ as a floating parameter *)
xsec1 = Subscript[\[Sigma], total] /. {gWWZ -> -0.7};
xsec2 = Subscript[\[Sigma], total] /. {gWWZ -> -0.9};
xsec3 = Subscript[\[Sigma], total] /. {gWWZ -> 0.9};
xsec4 = Subscript[\[Sigma], total] /. {gWWZ -> 0.7};

Plot[{xsec1, xsec2, xsec3, xsec4}, {Ecm, 160, 210}, Frame -> True, 
  PlotLegends -> {"\!\(\*SubscriptBox[\(g\), \(WWZ\)]\) = -0.7", 
    "\!\(\*SubscriptBox[\(g\), \(WWZ\)]\) = -0.9", 
    "\!\(\*SubscriptBox[\(g\), \(WWZ\)]\) = 0.9", 
    "\!\(\*SubscriptBox[\(g\), \(WWZ\)]\) = 0.7"}, PlotRange -> {0, 30}];
Labeled[%, {"\[Sigma] (pb)", 
  "\!\(\*SubscriptBox[\(E\), \(cm\)]\) (GeV)"}, {Left, Bottom}]
