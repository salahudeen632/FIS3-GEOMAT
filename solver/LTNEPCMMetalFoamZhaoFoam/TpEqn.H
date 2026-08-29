// PCM energy equation: Eq. (6), apparent-heat-capacity formulation.
surfaceScalarField rhoCpPhi
(
    "rhoCpPhi",
    fvc::interpolate(rhoP*CpP)*phi
);

fvScalarMatrix TpEqn
(
    fvm::ddt(epsilon*rhoP*CpP, Tp)
  + fvm::div(rhoCpPhi, Tp)
  - fvm::laplacian(kPe, Tp)
  + fvm::Sp(hps*aPS, Tp)
 ==
    hps*aPS*Ts
);

TpEqn.relax();
TpEqn.solve();
