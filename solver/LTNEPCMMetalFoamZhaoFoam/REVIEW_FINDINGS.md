# Review findings

## High impact corrections

1. `rhoCp` was a temporary field created inside `TpEqn.H`. Its previous-time
   state was therefore not reliable for `fvc::ddt(rhoCp)`.
2. The prior `Sp(ddt(rhoCp), Tp)` product correction used the current
   temperature in the implicit source. For Euler time stepping this produces
   the previous-time heat-capacity coefficient instead of the current
   apparent heat capacity.
3. `PrPCM` used the latent-heat-enhanced apparent heat capacity. This can raise
   Prandtl number into the thousands and artificially increase `hps`. The
   reviewed solver uses sensible heat capacity for `PrPCM`.
4. The product-rule thermal correction was numerically fragile and could drive
   `Tp` to large negative values. It has been replaced by a direct implicit
   Euler apparent-heat-capacity transient term.

## Model checks retained

- `kPe` is used in the PCM conduction equation.
- `kSe` is used in the foam conduction equation.
- `kP`, not `kPe`, is used in the Zukauskas correlation.
- The mushy source remains `mushyCu*(1-b)^2/(b^3+q)`.
- Buoyancy remains `-rho*beta*(Tp-Tref)*g`.

## Case information still required

Accurate validation also requires the complete case and log. Check:

- `deltaT = 3 K`
- identical temperatures and PCM properties for all metals
- only `rhoFoam`, `CpFoam`, and `kFoam` change
- no stale time directories or `postProcessing`
- `useViscosityDamping false` when the mushy source is active

## Dimensionless-time definition

The solver writes one requested definition to
`postProcessing/dimensionlessTime/dimensionlessTime.dat`:

- fixed `Tb`
- fixed liquid `CpPl` for `Ste`
- fixed liquid `alphaP = kPl/(rhoPl*CpPl)` for `Fo`

This matches the conventional melting Stefan number and avoids inserting the
latent-heat apparent-capacity peak into the plotting coordinate.
