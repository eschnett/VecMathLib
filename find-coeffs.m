(* -*-Mathematica-*- script to determine the polynomial coefficients
   to approximate functions
   
   2013-02-12 Erik Schnetter <eschnetter@perimeterinstitute.ca>
   
   Based on the "minimax" algorithm, but using least squares instead *)

funcs =
  {{Sin[2 Pi #]&, 0, 1/4, False, 1, 2, 10},
   {Cos[2 Pi #]&, 0, 1/4, False, 0, 2, 10},
   {Tan, 0, Pi/4, False, 1, 2, 20},
   (* {Cot, 0, Pi/2, False, 0, 2, 10}, *)
   {Log[(1+#)/(1-#)]&, 0, 1/3, False, 1, 2, 15},
   {Exp, 0, 1, False, 0, 1, 15},
   {ArcTan, 0, 1, False, 1, 2, 25}};

findcoeffs[func_, xmin_, xmax_, pade_, dmin_, dstep_, ndegrees_] :=
  Module[{prec, npts,
          degree,
          qs, q2x, x2q,
          funcq,
          funcqpade, polydenom,
          coeffs, powers, poly,
          norm,
          A, r, b,
          approx,
          error1, error2, error3, error},
         
         (* Working precision *)
         prec = 30;
         
         (* Number of test points *)
         npts = 1000; (*10000;*)
         
         degree = (dmin + ndegrees dstep) / If[pade, 2, 1];
         
         (* Test points in normalized interval [0,1] *)
         qs = Table[n/(npts-1), {n, 0, npts-1}];
         
         (* A (discrete) L2-norm based on the test points *)
         norm[f_] = Sqrt[Total[N[(Abs[f[#1]]^2 &) /@ qs, prec]] / Length[qs]];
         
         (* Transform q to x coordinate *)
         q2x[q_] = xmin + (xmax - xmin) q;
         x2q[x_] = (x - xmin) / (xmax - xmin);
         
         (* Function with rescaled input *)
         funcq[q_] = func[q2x[q]];
         
         polydenom[q_] =
         If[pade,
            (* Use denominator of Pade approximant as denominator of
               our approximant *)
            funcqpade[q_] = PadeApproximant[funcq[q], {q, 0, degree}];
            funcqpade[q][[1,1]],
            (* else *)
            (* Use a trivial denominator *)
            1];
         
         (* List of expansion coefficients *)
         coeffs = Table[c[i], {i, dmin, degree, dstep}];
         
         (* Corresponding list of powers of q *)
         powers[q_] = Table[If[i==0, 1, q^i], {i, dmin, degree, dstep}];
         
         (* Polynomial *)
         poly[q_] = coeffs . powers[q];
         
         (* We determine the expansion coefficients via a
            least-squares method *)
         A = N[Table[powers[q], {q, qs}], prec];
         r = N[(funcq[#1] polydenom[#1] &) /@ qs, prec];
         (* r = N[(Limit[funcq[q] polydenom[q], q->#1] &) /@ qs, prec]; *)
         b = LeastSquares[A, r];
         
         (* Define approximating polynomial using this solution *)
         approx[q_] = ((poly[q] /. MapThread[#1 -> #2 &, {coeffs, b}]) /
                       N[polydenom[q], prec]);
         
         (* Calculate three kinds of errors to check solution: *)
         (* (1) the (discrete) norm form above: *)
         error1 = norm[approx[#1] - funcq[#1] &];
         (* (2) A non-discrete L2-norm using an integral: *)
         error2 = Sqrt[NIntegrate[Abs[approx[q] - funcq[q]]^2, {q, 0, 1},
                                  WorkingPrecision -> prec]];
         (* (3) The maximum of the error: *)
         error3 = NMaxValue[{Abs[approx[q] - funcq[q]], q >= 0 && q <= 1}, {q},
                            WorkingPrecision -> prec];
         error = Max[error1, error2, error3];
         
         Print[{func, ndegrees, CForm[error]}];
         Write[outfile, {func, ndegrees, CForm[error],
                         CForm[HornerForm[approx[x2q[x]]]]}]];

findcoeffs2[func_, xmin_, xmax_, pade_, dmin_, dstep_, maxndegrees_] :=
  Do[findcoeffs[func, xmin, xmax, pade, dmin, dstep, ndegrees],
     {ndegrees, maxndegrees}];

(*
outfile = First[Streams["stdout"]];
(* outfile = OpenWrite[]; *)
findcoeffs[ArcTan, 0, 1, False, 1, 2, 25];
*)

outfile = OpenWrite["coeffs.out"];
Write[outfile, "(* Coefficients for function approximations *)"];
Map[findcoeffs2@@# &, funcs];
Close[outfile];
