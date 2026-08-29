// Update temperature-dependent PCM properties, LTNE conductivities and
// foam/PCM interfacial heat transfer.
{
    const scalar Tsol = Tm.value() - 0.5*phaseChangeDeltaT.value();
    const scalar Tliq = Tm.value() + 0.5*phaseChangeDeltaT.value();

    // Eq. (9): Gaussian Delta-Dirac distribution.  The paper applies it
    // only over [Tsol,Tliq].  Normalisation by erf(2) preserves exactly Lf
    // over this truncated interval because deltaT4 = phaseChangeDeltaT/4.
    const scalar deltaT4 = 0.25*phaseChangeDeltaT.value();
    const scalar truncatedGaussianIntegral = std::erf(2.0);
    const scalar invSqrtPiDeltaT4 =
        1.0
       /
        (
            std::sqrt(constant::mathematical::pi)
           *deltaT4
           *truncatedGaussianIntegral
        );

    const scalarField& TpI = Tp.primitiveField();
    const vectorField& UI = U.primitiveField();

    scalarField& alphaI = alpha.primitiveFieldRef();
    scalarField& deltaDiracI = deltaDirac.primitiveFieldRef();
    scalarField& rhoPI = rhoP.primitiveFieldRef();
    scalarField& CpPI = CpP.primitiveFieldRef();
    scalarField& CpSensibleI = CpSensible.primitiveFieldRef();
    scalarField& rhoCpI = rhoCp.primitiveFieldRef();
    scalarField& kPI = kP.primitiveFieldRef();
    scalarField& muPI = muP.primitiveFieldRef();
    scalarField& kPeI = kPe.primitiveFieldRef();
    scalarField& kSeI = kSe.primitiveFieldRef();
    scalarField& kEffI = kEff.primitiveFieldRef();
    scalarField& hpsI = hps.primitiveFieldRef();
    scalarField& ReI = RePCM.primitiveFieldRef();
    scalarField& PrI = PrPCM.primitiveFieldRef();
    scalarField& mushyCoeffI = mushyCoeff.primitiveFieldRef();

    forAll(TpI, celli)
    {
        const scalar T = TpI[celli];

        // Eq. (4): liquid fraction b(T).
        const scalar b = min
        (
            max((T - Tsol)/(Tliq - Tsol), scalar(0)),
            scalar(1)
        );
        alphaI[celli] = b;

        // Eq. (9): latent-heat distribution D(T).
        scalar D = 0.0;
        if (T >= Tsol && T <= Tliq)
        {
            const scalar x = (T - Tm.value())/deltaT4;
            D = invSqrtPiDeltaT4*std::exp(-sqr(x));
        }
        deltaDiracI[celli] = D;

        // Eqs. (7), (8) and (10).
        const scalar rho =
            rhoPs.value() + (rhoPl.value() - rhoPs.value())*b;

        const scalar cpSensible =
            CpPs.value()
          + (CpPl.value() - CpPs.value())*b;

        const scalar cp = cpSensible + Lf.value()*D;

        const scalar k =
            kPs.value() + (kPl.value() - kPs.value())*b;

        // Eqs. (2)-(5): Carman-Kozeny type mushy damping.
        const scalar Fmush = sqr(1.0 - b)/(pow3(b) + q);

        rhoPI[celli] = rho;
        CpPI[celli] = cp;
        CpSensibleI[celli] = cpSensible;
        rhoCpI[celli] = rho*cp;
        kPI[celli] = k;
        muPI[celli] = muPl.value()
            *(1.0 + (useViscosityDamping ? viscosityA*Fmush : 0.0));

        // Paper Eqs. (18)-(19), obtained from Boomsma-Poulikakos Eq. (17).
        kPeI[celli]  = boomsmaKeff(0.0,           k, epsilon);
        kSeI[celli]  = boomsmaKeff(kFoam.value(), 0.0, epsilon);
        kEffI[celli] = boomsmaKeff(kFoam.value(), k, epsilon);

        // Reynolds and Prandtl numbers used by Eq. (12).  U is the
        // volume-averaged/superficial velocity solved by Eqs. (2)-(3).
        const scalar ReLocal =
            rho*mag(UI[celli])*dp/max(muPI[celli], SMALL);

        const scalar PrLocal =
            cpSensible*muPI[celli]/max(k, SMALL);

        ReI[celli] = ReLocal;
        PrI[celli] = PrLocal;

        scalar hpsCorrelation = 0.0;

        if (hpsModel == "Zukauskas")
        {
            // Paper Eq. (12).
            if (ReLocal > SMALL && ReLocal <= 40.0)
            {
                hpsCorrelation =
                    0.76
                   *std::pow(ReLocal, 0.4)
                   *std::pow(PrLocal, 0.37)
                   *k/dp;
            }
            else if (ReLocal > 40.0 && ReLocal <= 1000.0)
            {
                hpsCorrelation =
                    0.52
                   *std::pow(ReLocal, 0.5)
                   *std::pow(PrLocal, 0.37)
                   *k/dp;
            }
            else if (ReLocal > 1000.0 && ReLocal <= 20000.0)
            {
                hpsCorrelation =
                    0.26
                   *std::pow(ReLocal, 0.6)
                   *std::pow(PrLocal, 0.37)
                   *k/dp;
            }
            else if (ReLocal > 20000.0)
            {
                hpsCorrelation =
                    0.26
                   *std::pow(20000.0, 0.6)
                   *std::pow(PrLocal, 0.37)
                   *k/dp;
            }
        }
        else if (hpsModel == "constant")
        {
            // Required for a no-flow LTNE calculation because Eq. (12)
            // gives hps = 0 when Re = 0.
            hpsCorrelation = hpsConstant.value();
        }

        // hpsMin is an optional numerical/stagnant coupling floor.  Set it
        // to zero for a strict Eq. (12) calculation.
        hpsI[celli] = max(hpsMin.value(), hpsCorrelation);

        mushyCoeffI[celli] = rho*mushyCu.value()*Fmush;
    }

    alpha.correctBoundaryConditions();
    deltaDirac.correctBoundaryConditions();
    rhoP.correctBoundaryConditions();
    CpP.correctBoundaryConditions();
    CpSensible.correctBoundaryConditions();
    rhoCp.correctBoundaryConditions();
    kP.correctBoundaryConditions();
    muP.correctBoundaryConditions();
    kPe.correctBoundaryConditions();
    kSe.correctBoundaryConditions();
    kEff.correctBoundaryConditions();
    hps.correctBoundaryConditions();
    RePCM.correctBoundaryConditions();
    PrPCM.correctBoundaryConditions();
    mushyCoeff.correctBoundaryConditions();
}
