(* ====================================================== *)
(* GitHub: https://github.com/yuriphestos/gradphestos-epp *)
(* Graduate Level Elementary Particle Physics             *)
(* Random Tasks for Fun                                   *)
(* ====================================================== *)
(* Task 2                                                 *)
(* Task 2.1: Z Boson Cross Section at the Tevatron        *)
(* ====================================================== *)


(* pdfs for up and down quarks *)
u[x_] := Piecewise[{{0, x < 0}, {1/0.38*x^(-0.52) (1 - x)^(4.11), 
    0 < x < 1}, {0, x > 1}}]
d[x_] := Piecewise[{{0, x < 0}, {1/0.79*x^(-0.65) (1 - x)^(5.08), 
    0 < x < 1}, {0, x > 1}}]


(* m_Z^2/Tevatron cm energy^2 in GeV *)
r := 91.18^2/1960^2
r

Plot[{u[x]*u[r/x]/x, d[x]*d[r/x]/x}, {x, 0, 1}]

(* integral over x entering p antip \[Rule] Z production *)
IntOverX1 = NIntegrate[u[x]*u[r/x]/x, {x, r, 1}]
IntOverX2 = NIntegrate[d[x]*d[r/x]/x, {x, r, 1}]

(* Z decay width to uu and dd in GeV *)
\[CapitalGamma]Zuu = 0.3;
\[CapitalGamma]Zdd = 0.4;


(* Prefactor *)
Prf = (4*\[Pi]^2)/(91.19*1960^2)

(* Total cross section of pp -> Z *)
CrxppZ1 = IntOverX1*\[CapitalGamma]Zuu;
CrxppZ2 = IntOverX2*\[CapitalGamma]Zdd;

CrxppZ = Prf*(CrxppZ1 + CrxppZ2)

(* Branching ratio of Z -> ee and \[Mu]\[Mu] *)
BrZll = (0.084*2)/2.494

(* Total cross section of pp -> ll in pb *)

Crxppll = CrxppZ*BrZll*0.3894*10^9

(* Dorigo Z->ee+Z->\[Mu]\[Mu] at Tevatron in pb *)
Dorigo = 267 + 246;

Crxppll/Dorigo
