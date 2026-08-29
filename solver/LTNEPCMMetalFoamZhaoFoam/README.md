# LTNEPCMMetalFoamPaperFormFoam v11

OpenFOAM Foundation v11 solver for a PCM-saturated open-cell metal foam using
local thermal non-equilibrium (LTNE):

- `Tp`: PCM temperature
- `Ts`: metal-foam temperature
- `U`, `p`: volume-averaged PCM flow variables when `solveFlow yes`

## Implemented equations

The solver follows the supplied main document:

- PCM equation: Eq. (6), using `kPe` and `+hps*aPS*(Ts-Tp)`.
- Foam equation: Eq. (11), using `kSe` and `-hps*aPS*(Ts-Tp)`.
- PCM density, apparent heat capacity and conductivity: Eqs. (7)-(10).
- Liquid fraction / mushy damping: Eqs. (4)-(5).
- Boomsma-Poulikakos conductivity model: Eqs. (13)-(19).
- Foam `K`, `CForch`, `dp`, and `aPS`: Eqs. (20)-(24).

`kEff` is written for diagnostics only. It is **not** used in both LTNE
energy equations.

## hps closure

The paper's Eq. (12) is implemented as:

```text
hpsModel Zukauskas;
```

This needs a nonzero flow velocity. If `solveFlow no`, then `RePCM=0` and Eq.
(12) gives `hps=0`, which decouples `Tp` and `Ts`.

For a no-flow LTNE calculation, use:

```text
hpsModel    constant;
hpsConstant 1000;  // W/(m2 K), user-calibrated value
hpsMin      0;
```

`hpsConstant` and `hpsMin` are entered as bare numerical values in W/(m2 K).
The solver applies:

```text
hps = max(hpsMin, hpsCorrelation)
```

## Additional corrections in this version

- Uses paper Eq. (6) in non-conservative material form. During each thermal
  correction, `rhoP*CpP` is frozen and only `Tp` is time differentiated:
  `epsilon*rhoP*CpP*dTp/dt + rhoP*CpP*(U dot grad(Tp))`.
- Uses `-rho*beta*(Tp-Tref)*g`, so hot PCM rises when `g` points downward.
- Reads `mushyCu` as `[1 -3 -1 0 0 0 0]` and does not multiply it by density.
- Uses a dimensionless Calmidi Forchheimer coefficient compatible with
  `rho*CForch*|U|/sqrt(K)`.
- Deliberately retains `q` directly in `b^3 + q`, as requested.
- `hpsMin` is applied to internal cells; it is no longer overwritten.
- Adds selectable `hpsModel` (`Zukauskas` or `constant`).
- Writes `RePCM` and `PrPCM` fields.
- Normalises the truncated Gaussian `D(T)` so that the latent-heat integral
  over the declared phase-change range is exactly `Lf`.
- Updates LTNE fixed-gradient thermal boundaries from `qBottom`, `hLoss`, and
  `Tinf` each thermal correction instead of relying on static loss gradients.
- Prints min/max values of `Tp`, `Ts`, `hps`, `RePCM`, and `PrPCM`.

## Boundary heat inputs

The case dictionary must provide:

```text
Tw      [0 0 0 1 0 0 0] 313.15;
hotWallThermalMode fixedTemperature;
qHotWall [1 0 -3 0 0 0 0] 0;
Tinf    [0 0 0 1 0 0 0] 298;
qBottom [1 0 -3 0 0 0 0] 1600;
hLoss   [1 0 -3 -1 0 0 0] 3;
```

`Tw` is the single prescribed hot-wall temperature. The solver assigns it to
both `Tp` and `Ts` on the patch selected by `hotWallPatch`, before every
thermal equation solve.

The selected patch must be declared as `fixedValue` in both `0/Tp` and
`0/Ts`. Their file values are only initial placeholders because the solver
overwrites them with `Tw`.

Do not select `bottomWall` as `hotWallPatch` when also using `qBottom`; one
patch cannot simultaneously use prescribed temperature and prescribed heat
flux.

For a prescribed heat flux on the selected hot wall, use:

```foam
hotWallPatch       leftWall;
hotWallThermalMode heatFlux;
qHotWall           [1 0 -3 0 0 0 0] 1800;
qBottom            [1 0 -3 0 0 0 0] 0;
hLoss              [1 0 -3 -1 0 0 0] 0;
```

Declare the selected patch as `fixedGradient` in both `0/Tp` and `0/Ts`:

```foam
leftWall
{
    type      fixedGradient;
    gradient  uniform 0;
}
```

The initial gradient is a placeholder. At runtime the solver sets
`gradient=qHotWall/(kPe+kSe)`, so the two LTNE conductive phase fluxes sum to
the requested total wall heat flux. Other adiabatic walls can remain
`zeroGradient`.

The solver applies a common `Tp`/`Ts` gradient so the LTNE phase fluxes sum to
the requested total boundary flux.

## Runtime Ste, Fo and tau output

The solver writes:

```text
postProcessing/dimensionlessTime/dimensionlessTime.dat
```

Add these optional entries to `constant/pcmMetalFoamProperties`:

```text
hotWallPatch          leftWall;
characteristicLength  0.04;     // metres
TbulkReference        293.15;   // kelvin; 0 uses the initial volume average
```

`Tw` is required and is listed in the boundary heat-input section above.

For comparison with the paper, the output uses fixed liquid PCM properties:

The `TbFixed` column is the fixed `TbulkReference` value used in
the Stefan-number calculation. The apparent heat capacity containing latent
heat is not used to construct the dimensionless plotting axis.

```text
CpSte = CpPl
alphaP = kPl/(rhoPl*CpPl)
Ste = CpSte*abs(Tw-TbFixed)/Lf
Fo = alphaP*time/H^2
tau = Ste*Fo
```

The package includes the digitized paper Fig. 4 benchmark:

```text
benchmark_fig4_digitized.csv
benchmark_fig4.png
```

Run the comparison from the case directory after copying the script and
benchmark CSV there, or provide their paths explicitly:

```bash
python3 plot_dimensionless_time.py --metal Copper
python3 plot_dimensionless_time.py --metal Aluminum
python3 plot_dimensionless_time.py --metal Nickel
```

Each plot contains only the selected paper benchmark and the numerical melt
fraction against solver-computed `tau`.

The buoyancy implementation follows the standard OpenFOAM gravity-vector
convention. With `g = (0 -9.81 0)`, warmer PCM receives an upward force.

## Inertial coefficient note

The supplied document prints Eq. (21) with a `dp^2` multiplier. That factor is
not used here because `CForch` must be dimensionless in the implemented source
`rho*CForch*|U|/sqrt(K)`.

## Build

```bash
cd /home/salu/OpenFOAM/salu-11/applications/solvers
cp -r LTNEPCMMetalFoamDoc3Foam .
cd LTNEPCMMetalFoamDoc3Foam
wclean
wmake
```

The executable is:

```text
$FOAM_USER_APPBIN/LTNEPCMMetalFoamPaperFormFoam
```
