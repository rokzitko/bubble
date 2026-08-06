(* Broyden mixing, Rok Zitko, rok.zitko@ijs.si, May 2009 *)
(* See D. D. Johnson, Phys. Rev. B 38, 12807 (1988). *)

(* -- OCCUPANCY CONTROL -- *)
(* __ CLEAN VERSION __ *)

Print["broyden.m, $Id: broyden.m,v 1.1 2009/05/05 07:42:23 rok Exp rok $"];

(* Length of vector *)
dim0 = Length[Import["dmft/1-Delta.dat"]];
dim = dim0 + 1; (* mu *)

col[v_] := Transpose[{v}];
row[v_] := {v};
scalar[v1_, v2_] := Flatten[v1] . Flatten[v2];
sparseIdentityMatrix[n_] := SparseArray[{{i_, i_} -> 1}, {n, n}];

Print["offset=", offset, " alpha=", alpha];

(* iter-offset = M *)
offset = Max[{offset, iter-M}];

Print["offset(M)=", offset];

w[0] = 0.01;
w[_] = 1;

(* --- *)

DF[m_] := DF[m] = (F[m + 1] - F[m])/Norm[Flatten[F[m + 1] - F[m]]];
Dn[m_] := Dn[m] = (mu[m + 1] - mu[m])/Norm[Flatten[F[m + 1] - F[m]]];

a[m_] := a[m] = Table[w[i] w[j] scalar[DF[i], DF[j]], {i, m-1}, {j, m-1}];
beta[m_] := beta[m] = Inverse[w[0]^2 IdentityMatrix[m-1] + a[m]];
beta[m_, k_, n_] := beta[m][[k, n]];

u[n_] := u[n] = alpha DF[n] + Dn[n];
c[k_, m_] := c[k, m] = scalar[DF[k], F[m]];

mu[m_] := mu[m] = mu[m - 1] + alpha F[m - 1] - 
 Sum[w[n] w[k] c[k, m-1] beta[m - 1, k, n] u[n], {k, 1, m - 2}, {n, 1, m - 2}];

(* --- *)

freqs = Import["dmft/" <> ToString[offset] <> "-Delta.dat"][[All, 1]];

(* Initial approximation *)

mu[m_ /; m < mres] := mu[m] = col @ Join[
  Import["dmft/" <> ToString[offset + m - 1] <> "-Delta.dat"][[All, 2]],
{ Import["dmft/" <> ToString[offset + m - 1] <> "-mu.dat"][[1,1]] }
];

prop[m_] := prop[m] = col @ Join[
  Import["Delta-" <> ToString[offset + m] <> ".dat"][[All, 2]],
{ Import["dmft/" <> ToString[offset + m - 1] <> "-occup.dat"][[1,1]] }
];

(* TARGET OCCUPANCY *)
n0 = goal;
deriv = 0.5;

(* System definition *)
F[m_] := F[m] = Module[{tmp},
  tmp = prop[m] - mu[m];
  tmp[[dim0+1]] = deriv ( prop[m][[dim0+1]] - n0 );
  tmp
];

mres = iter-offset+2;

Foccup = Flatten[F[mres-1]][[-1]];
Print["F occup = ", Foccup];

res = Norm[Flatten[F[mres-1]]];
Print["residual= ", res];

dim = Dimensions[a[mres-1]];
Print["dim=", dim];

sv = SingularValues[a[mres-1]] [[2]];
Print["singular values=", sv];

cond = sv[[1]]/sv[[-1]];
Print["condition number=", cond];

norm = Norm @ Flatten[a[mres-1]];
Print["norm=", norm];

{Foccup, res, dim, norm, cond} >>> "BROYDEN";

muexp[m_] := Transpose[{freqs, Take[ mu[m] // N // Flatten, {1, dim0}]}];

saveres[m_] := Export["Delta-" <> ToString[offset + m - 1] <> "-broyden.dat", muexp[m]];
savemu[m_] := Export["param.eps", {Take[mu[m] // N // Flatten, {dim0+1, dim0+1}]}, "Table"];

saveres[mres];
savemu[mres];

Print["DONE"];
