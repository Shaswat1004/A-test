(* Mathematica cross-check for the two-loop figure-eight vacuum-energy graph. *)
ClearAll[d, Nn, T, eps, Km, m, m2, mub, L, Ns, np, theta, n, r];

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

(* N_s-string generalization.  The relative-mode projector is
   P_AB = delta_AB - 1/N_s, so sum_A P_AA^2=(N_s-1)^2/N_s. *)
ProjectorFactor = FullSimplify[Ns ((Ns - 1)/Ns)^2];
ExpectedProjectorFactor = (Ns - 1)^2/Ns;
Print["N_s-string projector factor: ", ProjectorFactor];
Print["Projector factor matches: ",
  FullSimplify[ProjectorFactor - ExpectedProjectorFactor] === 0];

DeltaRhoNs = FullSimplify[
  np ProjectorFactor (np (d - 2) - 2 d) Km^2/(8 T d)
  ];
ExpectedNs = np (Ns - 1)^2 (np (d - 2) - 2 d) Km^2/
  (8 T Ns d);
Print["N_s-string general-d result: ", DeltaRhoNs];
Print["N_s-string general-d result matches: ",
  FullSimplify[DeltaRhoNs - ExpectedNs] === 0];

(* N_s=1 has no relative modes; N_s=2 must reproduce the old result. *)
Print["N_s=1 relative correction vanishes: ",
  FullSimplify[DeltaRhoNs /. Ns -> 1] === 0];
Print["N_s=2 reproduces the two-string result: ",
  FullSimplify[(DeltaRhoNs /. Ns -> 2) - (DeltaRho /. Nn -> np)] === 0];
Print["N_s=2 mass normalization: ",
  FullSimplify[-T m^2/(2 Ns) /. Ns -> 2] === -T m^2/4];

CoeffNs = Series[
   ExpectedNs/Km^2 /. d -> 2 - 2 eps,
   {eps, 0, 2}
   ] // Normal // Expand;
ExpectedCoeffNs = (
  -np (Ns - 1)^2/(4 T Ns)
  -np^2 (Ns - 1)^2 eps/(8 T Ns)
  -np^2 (Ns - 1)^2 eps^2/(8 T Ns)
  );
Print["N_s-string coefficient expansion: ", CoeffNs];
Print["N_s-string coefficient expansion matches: ",
  FullSimplify[CoeffNs - ExpectedCoeffNs] === 0];

BareSeriesNs = Series[
    (ExpectedNs/Km^2 /. d -> 2 - 2 eps) KmSeries^2,
    {eps, 0, 0}
    ] // Normal // Expand // FullSimplify;
ExpectedBareNs = -np (Ns - 1)^2 m^4/(64 Pi^2 T Ns) (
    1/eps^2 + (2 L + np/2)/eps
    + 2 L^2 + Pi^2/6 + np L + np/2
    );
Print["N_s-string bare expansion matches: ",
  FullSimplify[BareSeriesNs - ExpectedBareNs] === 0];

(* Explicit N_s=3 orthogonal basis and basis-independence check. *)
O3 = {
  {1/Sqrt[3], 1/Sqrt[2], 1/Sqrt[6]},
  {1/Sqrt[3], -1/Sqrt[2], 1/Sqrt[6]},
  {1/Sqrt[3], 0, -2/Sqrt[6]}
  };
Relative3 = O3[[All, 2 ;; 3]];
P3 = FullSimplify[Relative3 . Transpose[Relative3]];
ExpectedP3 = IdentityMatrix[3] - ConstantArray[1/3, {3, 3}];
Print["Three-string basis is orthogonal: ",
  FullSimplify[Transpose[O3] . O3 - IdentityMatrix[3]] ===
    ConstantArray[0, {3, 3}]];
Print["Three-string relative projector matches: ",
  FullSimplify[P3 - ExpectedP3] === ConstantArray[0, {3, 3}]];
Print["Three-string diagonal-square sum is 4/3: ",
  FullSimplify[Total[Diagonal[P3]^2]] === 4/3];

Rotation2 = {{Cos[theta], -Sin[theta]}, {Sin[theta], Cos[theta]}};
Relative3Rotated = Relative3 . Rotation2;
P3Rotated = FullSimplify[
  Relative3Rotated . Transpose[Relative3Rotated],
  Assumptions -> Element[theta, Reals]
  ];
Print["Relative-basis rotation leaves projector invariant: ",
  FullSimplify[P3Rotated - P3, Assumptions -> Element[theta, Reals]] ===
    ConstantArray[0, {3, 3}]];

(* General local binding potential.  Nn is the manuscript's component
   multiplicity N; the built-in symbol N is intentionally not used. *)
MomentGamma = 2^n Gamma[n + Nn/2]/Gamma[Nn/2];
MomentProduct = Product[Nn + 2 r, {r, 0, n - 1}];
MomentAssumptions = Element[n, Integers] && n >= 0 && Nn > 0;

Print["Potential moment gamma/product identity: ",
  FullSimplify[MomentGamma - MomentProduct,
    Assumptions -> MomentAssumptions] === 0];

LowPotentialMoments = Table[
  FullSimplify[MomentGamma /. n -> k, Assumptions -> Nn > 0],
  {k, 1, 4}
  ];
ExpectedLowPotentialMoments = {
  Nn,
  Nn (Nn + 2),
  Nn (Nn + 2) (Nn + 4),
  Nn (Nn + 2) (Nn + 4) (Nn + 6)
  };
Do[
  Print["Potential moment n=", k, ": ", LowPotentialMoments[[k]]],
  {k, 1, 4}
  ];
Print["Potential moments n=1,2,3,4 match: ",
  FullSimplify[LowPotentialMoments - ExpectedLowPotentialMoments] ===
    ConstantArray[0, 4]];

N1PotentialMoments = FullSimplify[LowPotentialMoments /. Nn -> 1];
Print["N=1 moments divided by I_m^n: ", N1PotentialMoments];
Print["N=1 moments give 1,3,15,105: ",
  N1PotentialMoments === {1, 3, 15, 105}];
Print["Potential-moment Wick recursion: ",
  FullSimplify[
    (MomentGamma /. n -> n + 1) - (Nn + 2 n) MomentGamma,
    Assumptions -> MomentAssumptions] === 0];
Print["Explicit q^4 Wick pairings give N(N+2): ",
  FullSimplify[Nn^2 + Nn + Nn - Nn (Nn + 2)] === 0];

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
