(* Mathematica cross-check for the two-loop figure-eight vacuum-energy graph. *)
ClearAll[d, Nn, T, eps, Km, m, m2, mub, L];

(* Coincident derivative propagator:
   <d_a phi^i d_b phi^j> = delta^{ij} delta_{ab} Km/d. *)

A2 = Nn (Nn d + 2) Km^2/d;
B2 = Nn (Nn + d + 1) Km^2/d;

DeltaRho = FullSimplify[A2/(16 T) - B2/(8 T)];
Expected = Nn (Nn (d - 2) - 2 d) Km^2/(16 T d);
Print["General-d Wick result: ", FullSimplify[DeltaRho]];
Print["Matches expected: ", FullSimplify[DeltaRho - Expected] === 0];

(* Strict d=2 coefficient. *)
Strict2D = FullSimplify[Expected /. d -> 2];
Print["Strict d=2 coefficient: ", Strict2D];
Print["Expected -N Km^2/(8T): ", FullSimplify[Strict2D + Nn Km^2/(8 T)] === 0];

(* N=1 check: in d=2, B=A^2 and L4,E = -A^2/(16T), while <A^2>=2 Km^2. *)
N1 = FullSimplify[Strict2D /. Nn -> 1];
Print["N=1 check: ", N1];

(* Keep d=2-2 eps until after multiplying the divergent integral. *)
CoeffEps = Series[
   Nn (Nn (d - 2) - 2 d)/(16 T d) /. d -> 2 - 2 eps,
   {eps, 0, 2}
   ] // Normal // Expand;
Print["Coefficient expansion: ", CoeffEps];

ExpectedCoeff = -Nn/(8 T) - Nn^2 eps/(16 T) - Nn^2 eps^2/(16 T);
Print["Coefficient expansion matches: ",
  FullSimplify[CoeffEps - ExpectedCoeff] === 0];

(* MSbar tadpole.  With measure
   (Exp[EulerGamma] mub^2/(4 Pi))^eps Integral[d^d p/(2Pi)^d],
   I_m = 1/(4Pi) Exp[gamma eps] Gamma[eps] (mub^2/m^2)^eps. *)
L = Log[mub^2/m^2];
ImSeries = 1/(4 Pi) (1/eps + L + eps (L^2/2 + Pi^2/12));
KmSeries = -m^2 ImSeries;

BareSeries = Series[
    (Nn (Nn (d - 2) - 2 d)/(16 T d) /. d -> 2 - 2 eps)
      KmSeries^2,
    {eps, 0, 0}
    ] // Normal // Expand // FullSimplify;

ExpectedBare = -Nn m^4/(128 Pi^2 T) (
    1/eps^2 + (2 L + Nn/2)/eps
    + 2 L^2 + Pi^2/6 + Nn L + Nn/2
    );

Print["Bare expansion through finite order: ", BareSeries];
Print["Bare expansion matches: ",
  FullSimplify[BareSeries - ExpectedBare] === 0];

(* Dimensional-analysis check in d=2: [m^4/T] = 4 - 2 = 2. *)
Print["Dimension check [m^4/T] = 2: ", 4 - 2 == 2];

(* One-loop determinant derivative check:
   d/d(m^2) [m^2/(8Pi) (1 + Log[mub^2/m^2])]
   = Log[mub^2/m^2]/(8Pi) = (1/2) I_m^ren. *)
oneLoopMSbar = m2/(8 Pi) (1 + Log[mub^2/m2]);
oneLoopDerivative = FullSimplify[D[oneLoopMSbar, m2]];
ExpectedHalfTadpole = Log[mub^2/m2]/(8 Pi);
Print["One-loop derivative: ", oneLoopDerivative];
Print["One-loop derivative matches half the renormalized tadpole: ",
  FullSimplify[oneLoopDerivative - ExpectedHalfTadpole] === 0];

(* Important conceptual checks, not algebraic identities:
   1. Mixed massless/massive vacuum contractions vanish only because the
      massless loop is scaleless in infinite-volume dimensional regularization.
   2. On a finite cylinder, the massless sum is not scaleless; do not drop it.
   3. The mass term breaks the relative shift symmetry of the NG action and
      must originate from additional physics; that physics may also generate
      further local interactions/counterterms.
   4. Do not quote the finite part of the bare figure-eight graph as a
      scheme-independent physical vacuum energy without specifying the full
      counterterm/renormalization prescription.
*)
