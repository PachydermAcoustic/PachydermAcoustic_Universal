# Experimental trihedral cabinet diffraction

Pach_Column_Array_Source offers Diffraction=TrihedralExperimental alongside Vanderkooy and NumericalHybrid. Recreate a column to regenerate its stored balloons. Universal exposes Cabinet_Diffraction_Trihedral with the existing Pressure and Driver_Balloon interfaces; its optional vertexScale=0 restores DED exactly. Universal uses Hare geometry only.

This experimental method adds four **front-corner** residuals to the existing DED field. It does not yet model all eight cabinet vertices. The original DED model and numerical hybrid remain available.

## Formulation

A rigid rectangular cabinet corner is locally the exterior of a solid octant, with acoustic Neumann conditions on three faces. Its angular domain on the unit sphere has three quarter-great-circle boundary arcs.

The canonical solver follows Bonner, Graham and Smyshlyaev, *The Computation of Conical Diffraction Coefficients in High-Frequency Acoustic Wave Scattering*, SIAM Journal on Numerical Analysis 43(3), 1202–1230 (2005), [DOI](https://doi.org/10.1137/040603358). Bonner's [2003 thesis](https://purehost.bath.ac.uk/ws/portalfiles/portal/188152548/Bradley_David_Bonner_thesis.pdf), Chapters 3–6, supplies the boundary-integral and spectral formulation. This is separate from Candy's whole-cabinet MFS method.

Classes_Cabinet_TrihedralReference.cs contains the spherical boundary discretization and imaginary-axis reference. Classes_Cabinet_TrihedralGeneral.cs supplies complex kernels, the general spectral contour and face traces. Classes_Cabinet_Trihedral.cs prepares the face-source residual and integrates it into DED.

The selected option uses four cubic-graded panels per arc, eight-point panel quadrature, spectral cutoff 16 and explicit Abel damping 0.75. Reciprocity moves the source-on-baffle direction to a continuous single-layer trace, integrated with squared-endpoint quadrature. A transpose solve prepares that trace once per incident angle. Regular parts of the kernel derivative are tabulated; analytic singular terms are retained. Complex angular residuals are sampled on a 25 × 48 grid and interpolated cubically.

The known damped flat-face image is subtracted before adding the residual to DED. For each front corner, the incident direction follows the actual driver-to-corner vector on the baffle. Source and outgoing spherical propagation retain complex phase. The reference's positive-exponential coefficient is conjugated to match Pachyderm's negative-exponential propagation.

These are **engineering assumptions**, rather than formulas derived from the canonical paper:

- Low-frequency activation is (kL)^4 / (16 + (kL)^4), with smooth size L = 1/sqrt(1/width² + 1/height² + 1/depth²).
- Driver-to-corner distance receives an additional (kr)^2 / (16 + (kr)^2) activation.
- Driving pressure uses the existing grazing circular-piston DED factor.
- A cubic exterior taper spans direction cosine 0 through 0.05. Angular samples inside or on the closed solid octant are zero. Interpolation and the taper supply a continuous experimental extension near the faces.

These choices enable testing; they are not a uniform physical edge/vertex transition solution.

## Checks

Run from the repository root:

    dotnet run --project diagnostics/CabinetDiffraction -- --trihedral <output-directory>
    dotnet run --project diagnostics/CabinetDiffraction -- --trihedral-column
    dotnet run --project diagnostics/CabinetDiffraction -- --trihedral-general
    dotnet run --project diagnostics/CabinetDiffraction -- --trihedral-cabinet <output-directory>

Restricted canonical checks include an independent [Mehler–Dirichlet integral, DLMF 14.12.1](https://dlmf.nist.gov/14.12.E1), kernel derivatives, manufactured harmonic fields, boundary refinement, spectral refinement, reciprocity and axis permutations. The imaginary-axis coefficient is restricted to a conservative M1 convergence certificate; that is not a cabinet visibility switch.

Measured reference fixtures:

- Independent kernel maximum relative difference: 6.83E-11.
- Worst 64-panel-per-arc Neumann manufactured-field error at tau=8: 0.184%.
- Neumann axial/15-degree coefficient: +0.0241981434i. Doubling its 32-panel mesh changes it by 5.68E-6.
- Six Dirichlet coefficients agree with Bonner Table 6.7 within 0.66%. This is a Dirichlet benchmark, not an independent published Neumann comparison.
- Damped general-contour versus imaginary-axis agreement in M1: 1.70E-6. Extending cutoff 16 to 24 in the tested M2 direction changes the result by 2.65E-6.
- Exact source-on-face trace: -0.0990744436 + 0.0026565490i. Exterior offsets 0.01, 0.001 and 0.0001 approach it with differences 8.94E-4, 8.89E-5 and 8.56E-6.
- Prepared face residual versus independent reciprocity trace after subtracting the flat-face image: 3.01E-10.

Cabinet checks cover finite full-angle horizontal pressure in all eight bands, centered symmetry, grazing continuity, exact zero-scale DED recovery, a nonzero selected correction and finite eight-band balloon strings. CSVs record complex angular fields and side/oblique/rear depth corrections. Depth smoothness is checked by halving the sweep spacing: maximum complex correction second difference falls from 2.25E-4 to 6.78E-5 (ratio 0.301). The largest unit-pressure correction in the horizontal fixture is 0.0653. Visual inspection of those cuts showed no new isolated spikes. These fixtures do not establish physical accuracy for arbitrary cabinet dimensions or exclude every possible angular artifact.

## Remaining limitations

Finite damping smooths the spectral field; the zero-damping physical limit is unverified. General-direction boundary-mesh refinement, angular-table refinement and near-field accuracy need further testing. The general reference rejects double-boundary traces; the cabinet uses the explicit engineering extension above in that region.

Subtracting the flat-face image does not remove all edge endpoint content already present in finite DED integrals. Uniform endpoint matching remains unfinished, so double-counting and amplitude errors are possible. Rear-vertex driving, incident directions from higher-order edge paths, all-eight-corner closure, ports and mechanical driver loading are future work. Depth dependence in the added term currently comes from smooth activation; it does not represent rear-corner propagation.

Stored balloons retain magnitude-only per-octave normalization. Use Pressure or diagnostic CSVs for complex absolute comparisons. Rhino menu selection and document creation still need interactive testing after loading the rebuilt plugin.

For the existing DED lineage, see Vanderkooy, *A Simple Theory of Cabinet Edge Diffraction*, JAES 39(12), 923–933 (1991), [AES](https://aes2.org/publications/elibrary-page/?id=5952), and Urban et al., *The Distributed Edge Dipole Model for Cabinet Diffraction Effects*, JAES 52(10), 1043–1059 (2004), [AES](https://aes.org/publications/elibrary-page/?id=13024).

## Built test plugins

Both Universal and Rhino targets compile with zero errors (existing repository warnings remain). The isolated builds keep dependencies beside the plugin:

- .NET 7: C:/Users/Arthu/.codex/tmp/cabinet-trihedral/plugin-net7/Pachyderm_Acoustic.rhp
- .NET Framework: C:/Users/Arthu/.codex/tmp/cabinet-trihedral/plugin-net48/Pachyderm_Acoustic.rhp

Load the appropriate rebuilt plugin with its companion DLLs, then run Pach_Column_Array_Source and choose Diffraction=TrihedralExperimental. The running Rhino session has not been changed automatically. Existing columns need recreation to use the selected method.

## Default-column quadrature regression

The command defaults (8 drivers, spacing 0.075 m, diameter 0.05 m, width 0.08 m, height 0.575 m, depth 0.10 m) exposed a face-trace roundoff bug on driver 3. Distinct quadrature points could have a dot product rounded to exactly one. The logarithmic kernel then rejected them as coincident.

Face-trace quadrature now supplies sin²(angle/2) from the angular difference directly. The general off-face trace includes the normal separation as well. This retains the small positive separation without losing it in subtraction from a rounded dot product.

The --trihedral-column regression checks logarithmic values and derivatives below dot-product resolution, rejection of actual coincidence, and all eight drivers' complete finite eight-band balloons. It passes. The existing general-contour and cabinet continuity/depth checks also pass. This fixes evaluation of valid default geometry; the experimental physical assumptions above are unchanged.
