// Volume-averaged Brinkman-Forchheimer momentum equation, Eqs. (2)-(3).
surfaceScalarField rhoPhi
(
    "rhoPhi",
    fvc::interpolate(rhoP/sqr(epsilon))*phi
);

volScalarField porousResistance
(
    "porousResistance",
    muP/K + rhoP*CForch*mag(U)/sqrt(K) + mushyCoeff
);

volVectorField buoyancy
(
    "buoyancy",
    rhoP*beta*(Tp - Tref)*g
);

fvVectorMatrix UEqn
(
    fvm::ddt(rhoP/epsilon, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(muP/epsilon, U)
  + fvm::Sp(porousResistance, U)
);

UEqn.relax();

if (piso.momentumPredictor())
{
    solve
    (
        UEqn
     ==
       -fvc::grad(p_rgh)
      + buoyancy
    );
}
