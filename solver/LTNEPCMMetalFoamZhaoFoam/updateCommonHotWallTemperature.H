{
    if (hotWallThermalMode != "fixedTemperature")
    {
        // Heat-flux mode is updated in updateThermalBoundaryGradients.H.
        // Mode "none" leaves all temperature BCs exactly as supplied.
    }
    else
    {
    if
    (
        !isA<fixedValueFvPatchScalarField>
        (
            Tp.boundaryField()[hotWallPatchi]
        )
     || !isA<fixedValueFvPatchScalarField>
        (
            Ts.boundaryField()[hotWallPatchi]
        )
    )
    {
        FatalErrorInFunction
            << "Patch " << hotWallPatchName
            << " must use fixedValue for both Tp and Ts."
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField& TpHotWall =
        refCast<fixedValueFvPatchScalarField>
        (
            Tp.boundaryFieldRef()[hotWallPatchi]
        );

    fixedValueFvPatchScalarField& TsHotWall =
        refCast<fixedValueFvPatchScalarField>
        (
            Ts.boundaryFieldRef()[hotWallPatchi]
        );

    TpHotWall == Tw.value();
    TsHotWall == Tw.value();
    TpHotWall.evaluate();
    TsHotWall.evaluate();
    }
}
