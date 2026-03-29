(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 4                                                 *)
(* Task 4.4: Majorana Neutrinos and e- e- -> W- W-        *)
(* ====================================================== *)

<< FeynCalc`

(* ================================================================== *)
(* Matrix Elements for t- and u-channels                              *)
(* ================================================================== *)
(* 
   Amplitude structure (only m_\[Nu] piece survives chirality argument):
   Mt = (-gW^2 m\[Nu] / 8t) \[Epsilon]1^\[Mu] \[Epsilon]2^\[Lambda] 
        [vbar \[Gamma]\[Mu] (1+\[Gamma]5) \[Gamma]\[Lambda] (1-\[Gamma]5) u]
   Mu = same with k1 <-> k2 (\[Epsilon]1 <-> \[Epsilon]2)
*)

Mt = PolarizationVector[k1, \[Mu]] . PolarizationVector[k2, \[Lambda]] . 
     SpinorVBar[p1, 0] . GA[\[Mu]] . (1 + GA5) . GA[\[Lambda]] . (1 - GA5) . 
     SpinorU[p2, 0]

Mu = PolarizationVector[k2, \[Mu]] . PolarizationVector[k1, \[Lambda]] . 
     SpinorVBar[p1, 0] . GA[\[Mu]] . (1 + GA5) . GA[\[Lambda]] . (1 - GA5) . 
     SpinorU[p2, 0]

(* ================================================================== *)
(* Spin-summed, polarization-summed |M|^2 for each piece              *)
(* ================================================================== *)

(* |Mt|^2 *)
Mt ComplexConjugate[Mt];
FermionSpinSum[%];
DoPolarizationSums[%, k1];
DoPolarizationSums[%, k2];
Mtsq = (gW^4 Subscript[m, \[Nu]]^2 / (64 t^2)) * 
       Simplify[DiracSimplify[%]]

(* |Mu|^2 *)
Mu ComplexConjugate[Mu];
FermionSpinSum[%];
DoPolarizationSums[%, k1];
DoPolarizationSums[%, k2];
Musq = (gW^4 Subscript[m, \[Nu]]^2 / (64 u^2)) * 
       Simplify[DiracSimplify[%]]

(* Interference: Mt* Mu + Mu* Mt = 2 Re[Mt* Mu] *)
ComplexConjugate[Mt] Mu;
FermionSpinSum[%];
DoPolarizationSums[%, k2];
DoPolarizationSums[%, k1];
Mtusq = (gW^4 Subscript[m, \[Nu]]^2 / (64 t u)) * 
        Simplify[DiracSimplify[%]]

(* ================================================================== *)
(* Substitute Mandelstam variables and kinematics                     *)
(* ================================================================== *)
(* 
   On-shell conditions: p1^2 = p2^2 = 0 (massless electrons),
                        k1^2 = k2^2 = mW^2
   Mandelstam: s + t + u = 2 mW^2 (for me = 0)
   CM frame:   t = mW^2 - (s/2)(1 - \[Beta] Cos\[Theta])
               where \[Beta] = Sqrt[1 - 4mW^2/s]
*)

ScalarProduct[p1, p1] = 0;
ScalarProduct[p2, p2] = 0;
ScalarProduct[k1, k1] = mW^2;
ScalarProduct[k2, k2] = mW^2;
ScalarProduct[p1, p2] = s/2;
ScalarProduct[k1, k2] = s/2;  (* from (k1+k2)^2 = s and k1^2=k2^2=mW^2 *)
ScalarProduct[p1, k1] = (mW^2 - t)/2;
ScalarProduct[p2, k2] = (mW^2 - t)/2;
ScalarProduct[p1, k2] = (mW^2 - u)/2;
ScalarProduct[k1, p2] = (mW^2 - u)/2;

Eps[Momentum[k1], Momentum[k2], Momentum[p1], Momentum[p2]] = 0;

s = Ecm^2;
u = 2*mW^2 - s - t;
t = mW^2 - (s/2)*(1 - Sqrt[1 - (4*mW^2/s)]*Cos[\[Theta]]);

(* Total |M|^2 summed over spins and polarizations *)
Msq = Simplify[Mtsq + Musq + 2*Mtusq]

(* ================================================================== *)
(* Differential and Total Cross Section                               *)
(* ================================================================== *)
(* 
   d\[Sigma]/d\[CapitalOmega] = (1/(64\[Pi]^2 s)) * (|k1|/|p1|) * |M|^2
   |k1|/|p1| = Sqrt[1 - 4mW^2/Ecm^2]
*)

\[Sigma] = Simplify[
  Integrate[
    (1/(64*\[Pi]^2*Ecm^2)) * 
    Sqrt[1 - 4*mW^2/Ecm^2] * Msq * Sin[\[Theta]],
    {\[Theta], 0, \[Pi]},
    {\[CurlyPhi], 0, 2*\[Pi]},
    Assumptions -> {Ecm > 0 && mW > 0 && Ecm > 2*mW}
  ]
]

Print[""]
Print["==================================================="]
Print["Total cross section (symbolic, in terms of gW, m\[Nu], mW, Ecm):"]
Print[\[Sigma]]
Print["==================================================="]

(* ================================================================== *)
(* Numerical Evaluation                                               *)
(* ================================================================== *)
(* 
   gW = 2 mW / v where v = 246.22 GeV (Higgs vev)
   Equivalently: gW^4 / (64 mW^4) = GF^2 / 2
   We use gW = 2*80.377/246.22 \[TildeTilde] 0.6529
*)

Print[""]
Print["==================================================="]
Print["NUMERICAL EVALUATION"]
Print["==================================================="]

gW = 2*80.377/246.22;
Print["gW = ", gW]

(* --- At Ecm = 500 GeV (ILC-scale e-e- collider) --- *)
mW = 80.377;
Subscript[m, \[Nu]] = 0.8*10^(-9);  (* 0.8 eV in GeV *)
Ecm = 500;

\[Sigma]500 = \[Sigma];
Print[""]
Print["--- Ecm = 500 GeV ---"]
Print["\[Sigma] = ", \[Sigma]500, " GeV^(-2)"]
Print["\[Sigma] = ", \[Sigma]500*0.3894*10^12, " fb"]
Print["\[Sigma] = ", \[Sigma]500*0.3894*10^(-27), " cm^2"]
Print["N_events (L = 10^4 fb^-1) = ", 10^4 * \[Sigma]500*0.3894*10^12]

(* --- At Ecm = 13000 GeV (for comparison) --- *)
ClearAll[Ecm]
Ecm = 13000;

\[Sigma]13 = \[Sigma];
Print[""]
Print["--- Ecm = 13 TeV ---"]
Print["\[Sigma] = ", \[Sigma]13, " GeV^(-2)"]
Print["\[Sigma] = ", \[Sigma]13*0.3894*10^12, " fb"]
Print["\[Sigma] = ", \[Sigma]13*0.3894*10^(-27), " cm^2"]
Print["N_events (L = 140 fb^-1) = ", 140 * \[Sigma]13*0.3894*10^12]
