(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 1                                                 *)
(* Task 1.3: Near-Resonance Cross Section                 *)
(* ====================================================== *)


(* Constants in GeV -- Using standard PDG/Tree-level values *)
mZ = 91.1876;
\[CapitalGamma]Zhad = 1.744;
\[CapitalGamma]Zee = 0.0839;
\[CapitalGamma]Znu = 0.166;
\[Alpha] = 1/137.036;
me = 0.000511;

(* Conversion factor: 1/GeV^2 to nb *)
GeV2tonb = 3.89379*10^5;

(* Total Z width as a function of neutrino generations *)
\[CapitalGamma]Z[Nnu_] := \[CapitalGamma]Zhad + 3*\[CapitalGamma]Zee +
    Nnu*\[CapitalGamma]Znu;

(* Relativistic Breit-Wigner with s-dependent width *)
Subscript[\[Sigma], had][s_, 
   Nnu_] := (12*Pi*
     s/mZ^2)*(\[CapitalGamma]Zee*\[CapitalGamma]Zhad/((s - 
           mZ^2)^2 + (s^2/mZ^2)*\[CapitalGamma]Z[Nnu]^2))*GeV2tonb;

(* Full Kuraev-Fadin ISR Radiator Function Convolution *)
(* Includes the (1+3\[Beta]/4) virtual/soft photon correction factor *)
Subscript[\[Sigma], ISR][Ecm_?NumericQ, Nnu_] := 
 Module[{s = 
    Ecm^2, \[Beta]}, \[Beta] = (2*\[Alpha]/Pi)*(Log[s/me^2] - 1);
  NIntegrate[((1 + 3*\[Beta]/4) - 
     y^(1/\[Beta] - 1) + (1/2)*y^(2/\[Beta] - 1))*
   Subscript[\[Sigma], had][s*(1 - y^(1/\[Beta])), Nnu], {y, 0, 1}, 
  Method -> "GlobalAdaptive", PrecisionGoal -> 4]]

Plot[{Subscript[\[Sigma], ISR][Ecm, 2], 
  Subscript[\[Sigma], ISR][Ecm, 3], 
  Subscript[\[Sigma], ISR][Ecm, 4]}, {Ecm, 85, 95}, 
 PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Red], 
   Directive[Thick, Darker[Green]]}, 
 PlotLegends -> {"\!\(\*SubscriptBox[\(N\), \(\[Nu]\)]\) = 2", 
   "\!\(\*SubscriptBox[\(N\), \(\[Nu]\)]\) = 3", 
   "\!\(\*SubscriptBox[\(N\), \(\[Nu]\)]\) = 4"}, PlotRange -> All, 
 PlotTheme -> "Detailed", GridLines -> Automatic, Frame -> True, 
 FrameLabel -> {"\!\(\*SubscriptBox[\(E\), \(cm\)]\) (GeV)", 
   Rotate["\!\(\*SubscriptBox[\(\[Sigma]\), \(had\)]\) (nb)", 0]}]
