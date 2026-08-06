(* New approach: Sigma contains the ed shift! *)
(* Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Oct 2008 *)

(* Tabulated self-energy function *)
sigma = Import["res/resigma.dat"];
sigim = Import["res/imsigma.dat"];
sigma[[All, 2]] += I * sigim[[All, 2]];

fsigma = Interpolation[sigma];

(* Hilbert transform of the DOS *)
htDOS0[z_] := 2  (z - I Sign[Im[z]] Sqrt[1 - z^2]);
EPS = 10^-20;
htDOS[z_] := htDOS0[Re[z] + I (If[Im[z] > 0.0, Im[z], EPS])];

gf = Map[{#[[1]], htDOS[ #[[1]] - #[[2]] ]}&, sigma];
fgf = Interpolation[gf];

gfim = Map[{#[[1]], Im[#[[2]]]}&, gf];
gfre = Map[{#[[1]], Re[#[[2]]]}&, gf];

awim = Map[{#[[1]], -1/Pi #[[2]] }&, gfim];
awre = Map[{#[[1]], -1/Pi #[[2]] }&, gfre];

Export["res/imaw.dat", awim];
Export["res/reaw.dat", awre];
