// PISO pressure correction for p_rgh in Pa.
while (piso.correct())
{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA
    (
        constrainHbyA(rAU*(UEqn.H() + buoyancy), U, p_rgh)
    );

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    adjustPhi(phiHbyA, U, p_rgh);
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    while (piso.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, pRefValue);
        p_rghEqn.solve();

        if (piso.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();
        }
    }

    #include "continuityErrs.H"

    U = HbyA - rAU*fvc::grad(p_rgh);
    U.correctBoundaryConditions();

    p = p_rgh + rhoP*gh;
    p.correctBoundaryConditions();
}
