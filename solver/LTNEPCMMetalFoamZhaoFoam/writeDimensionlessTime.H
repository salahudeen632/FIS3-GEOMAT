{
    const scalarField& cellVolume = mesh.V();

    const scalar meltFraction =
        gSum(cellVolume*alpha.primitiveField())/totalVolume;

    scalar TwWall = Tw.value();

    if (hotWallPatchi >= 0)
    {
        const scalarField& hotWallArea =
            mesh.boundary()[hotWallPatchi].magSf();
        const scalarField& hotWallTemperature =
            Tp.boundaryField()[hotWallPatchi];

        TwWall = gSum(hotWallArea*hotWallTemperature)/gSum(hotWallArea);
    }
    else
    {
        TwWall = -GREAT;
        forAll(Tp.boundaryField(), patchi)
        {
            TwWall = max(TwWall, gMax(Tp.boundaryField()[patchi]));
        }
    }

    const scalar CpSte = CpPl.value();

    const scalar alphaP =
        kPl.value()/(rhoPl.value()*CpPl.value());

    const scalar Ste =
        CpSte*mag(TwWall - fixedBulkTemperature)/Lf.value();

    const scalar Fo =
        alphaP*runTime.value()/sqr(characteristicLength);

    const scalar tau = Ste*Fo;

    if (Pstream::master())
    {
        dimensionlessTimeFilePtr()
            << runTime.value() << ' '
            << fixedBulkTemperature << ' '
            << TwWall << ' '
            << CpSte << ' '
            << alphaP << ' '
            << meltFraction << ' '
            << Ste << ' '
            << Fo << ' '
            << tau << endl;
    }
}
