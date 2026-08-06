def1ch[1];

If[!paramexists["spin", "extra"],
   MyError["Define the spin of the impurity!"];
];
SPIN = ToExpression @ param["spin", "extra"];
MyPrint["SPIN=", SPIN];

Module[{sz, sp, sm, sx, sy, oz, op, om, ss},
       sz = spinketbraZ[SPIN];
       sp = spinketbraP[SPIN];
       sm = spinketbraM[SPIN];
       sx = spinketbraX[SPIN];
       sy = spinketbraY[SPIN];

       oz = nc[ sz, spinz[ d[] ] ];
       op = nc[ sp, spinminus[ d[] ] ];
       om = nc[ sm, spinplus[ d[] ] ];

       ss = oz + 1/2 (op + om) // Expand;

       H1 = eps number[d[]] + U hubbard[d[]] + B spinz[d[]];
       Hk = Jkondo ss + B sz; (* B on f-site *)

       H = H0 + H1 + Hc + Hk;
       Himp = H1 + Hk; (* Includes field and eps *)
       Hpot = U hubbard[d[]] + Jkondo ss; (* No field! *)
];

MAKESPINKET = SPIN;

(* All operators which contain d[], except hybridization (Hc). *)
Hselfd = H1 + Hk;

selfopd = ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];
