# Codex task: verify and maintain the vacuum-energy correction

Work only on the vacuum-energy calculation associated with Eq. (11) / `eq:two-string-ng-sum-relative` in this repository. Do not change unrelated physics or notation.

## Physics assumptions to preserve

1. The starting point is the two-string Nambu--Goto Lagrangian in `main.tex`, especially `eq:two-string-ng-sum-relative`.
2. The sum modes `X_s^i` are massless.
3. The relative modes `X_r^i` are assigned a common mass `m` by adding `-(T/4) m^2 X_r^i X_r^i` in Minkowski signature.
4. This mass deformation is external to Eq. (11), because Eq. (11) has a relative shift symmetry. State this caveat explicitly.
5. Let `N` be the number of massive relative components; only set `N=D-2` if every transverse relative mode is massive.
6. Use dimensional regularization in `d=2-2 epsilon` and infinite flat worldsheet volume for the figure-eight calculation.
7. Do not silently extend the infinite-volume dimensional-regularization result to a finite cylinder. On a cylinder, massless sums are not scaleless and mixed terms generally survive.

## Derivation that must be independently checked

Canonical normalization is

`phi_s^i = sqrt(T/2) X_s^i`, `phi_r^i = sqrt(T/2) X_r^i`.

The quadratic Minkowski theory is

`L2 = -1/2 (partial phi_s)^2 -1/2 (partial phi_r)^2 -1/2 m^2 phi_r^2`.

After Wick rotation, the pure massive quartic interaction should be

`L4E,r = +(1/(16 T)) (partial_a phi^i partial_a phi^i)^2
         -(1/(8 T)) (partial_a phi^i partial_b phi^i)
                     (partial_a phi^j partial_b phi^j)`.

Check the Wick-rotation signs from first principles rather than trusting this line.

Define

`K_m = Integral_p p^2/(p^2+m^2)`

with the MSbar measure. Rotational/O(N) invariance implies

`<partial_a phi^i partial_b phi^j> = delta^ij delta_ab K_m/d`.

Let

`A = partial_a phi^i partial_a phi^i`

and

`B = (partial_a phi^i partial_b phi^i)(partial_a phi^j partial_b phi^j)`.

Enumerate all three Wick pairings explicitly. The expected results are

`<A^2> = N(N d + 2) K_m^2/d`

and

`<B> = N(N+d+1) K_m^2/d`.

Therefore the connected two-loop one-vertex vacuum graph should be

`Delta rho_fig8 = <L4E,r>
 = N [N(d-2)-2d] K_m^2/(16 T d)`.

At strict d=2 the coefficient is `-N/(8T)`, but DO NOT set `d=2` before multiplying by the divergent integral. The `N(d-2)` term is evanescent and changes the single-pole and finite terms.

Use

`K_m = -m^2 I_m`,

`I_m = Integral_p 1/(p^2+m^2)`.

With the convention

`Integral_p = (Exp[gamma_E] mu_bar^2/(4 Pi))^epsilon Integral d^d p/(2 Pi)^d`,

check

`I_m = 1/(4 Pi) [1/epsilon + L + epsilon(L^2/2 + Pi^2/12) + O(epsilon^2)]`,

where `L=log(mu_bar^2/m^2)`.

The expected bare graph through finite order is

`Delta rho_fig8^bare = - N m^4/(128 Pi^2 T) [
  1/epsilon^2
  + (2 L + N/2)/epsilon
  + 2 L^2 + Pi^2/6 + N L + N/2
  + O(epsilon)
]`.

Do not call the finite term above a physical prediction without specifying the complete counterterm prescription. A local vacuum-energy/string-tension counterterm is allowed, and the phenomenological mass deformation may require parameter-counterterm insertions into the one-loop determinant.

## One-loop ordering check

Distinguish the Gaussian one-loop determinant from the first interaction correction:

`rho_1loop = (N/2) Integral_p log(p^2+m^2)`.

After MSbar subtraction, up to an arbitrary vacuum-energy counterterm,

`rho_1loop^MSbar = N m^2/(8 Pi) [1 + log(mu_bar^2/m^2)]`.

The quartic Nambu--Goto interaction first contributes through the two-loop figure-eight graph at order `1/T`.

## Required independent checks

Perform all of the following before editing the manuscript:

1. **Original-field check**: set `X_s=0`, hence `X_1=X_r/2`, `X_2=-X_r/2`, and derive the pure-relative quartic terms directly from the two original Nambu--Goto copies. Confirm the factors `-T/64` and `+T/32` before canonical normalization.
2. **N=1 check**: for one scalar, verify `B=A^2`. In d=2 this gives `L4E=-A^2/(16T)` and `<A^2>=2 K_m^2`, hence `Delta rho=-K_m^2/(8T)`.
3. **Index-combinatorics check**: independently enumerate the three Wick contractions for `A^2` and `B`; do not infer them only from the final formula.
4. **Mass dimensions**: in two dimensions `[T]=2`, `[m]=1`, so `m^4/T` has energy-density dimension two.
5. **Large-T limit**: interaction correction must vanish as `T -> infinity`.
6. **Massless dim-reg limit**: the infinite-volume figure-eight graph vanishes as `m -> 0` because all tadpoles become scaleless.
7. **Finite-volume warning**: explicitly state that the previous item is not a statement about finite-size/Casimir energy.
8. **Evanescent-term check**: expand the general-d coefficient through at least `epsilon^2`; verify
   `-N/(8T) - N^2 epsilon/(16T) - N^2 epsilon^2/(16T) + ...`.
9. **One-loop derivative check**: differentiate the claimed MSbar one-loop expression with respect to `m^2`; it should equal one half of the renormalized tadpole.
10. **Sign check**: trace the relation between Minkowski `L_M`, Euclidean `L_E`, and the vacuum energy `W/V = <L4E>` at first order. Make sure no extra minus sign from the vacuum functional is missed.

## Mathematica verification

Run or inspect `vacuum_energy_crosscheck.wl`. If Mathematica is installed, execute it with `wolframscript -file vacuum_energy_crosscheck.wl`. Every Boolean comparison printed by the script should be `True`.

If Mathematica is not installed, reproduce the same symbolic checks in Python/SymPy or another CAS and record that as a substitute; do not claim Mathematica was run when it was not.

## LaTeX validation

The derivation lives in `vacuum_energy_correction.tex` and is included from `main.tex` immediately before `\\end{document}`.

Run, if available:

`latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`

If `latexmk` is unavailable, run `pdflatex -interaction=nonstopmode -halt-on-error main.tex` at least twice. Resolve actual LaTeX errors introduced by this change. Do not spend time fixing unrelated pre-existing warnings unless they prevent compilation.

Also grep the log for `Undefined control sequence`, `LaTeX Error`, `Emergency stop`, and duplicate labels.

## Editing constraints

- Preserve the existing notation and RevTeX structure.
- Keep the detailed derivation in `vacuum_energy_correction.tex`; do not bloat unrelated sections.
- Do not change Eq. (11) itself unless an independent algebraic check proves it wrong.
- Do not claim the mass term is generated by Nambu--Goto dynamics.
- Do not replace the regulated general-d result by the strict d=2 expression before renormalization.
- Do not present the bare finite remainder as scheme independent.
- If a check fails, stop and explain the discrepancy instead of forcing the expected answer.

## Git workflow

Work on the existing branch `agent/vacuum-energy-correction` if present. Review `git diff main...HEAD` before committing. Commit only files relevant to this calculation. Push the branch and update/open a draft PR against `main`. In the PR description list: physics assumptions, derived result, regulator dependence, checks performed, and any checks that could not be run locally.
