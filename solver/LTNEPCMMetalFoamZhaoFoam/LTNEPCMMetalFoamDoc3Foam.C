/*---------------------------------------------------------------------------*\
  LTNEPCMMetalFoamZhaoFoam
  PCM saturated in open-cell metal foam.
  OpenFOAM Foundation v11.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "uniformDimensionedFields.H"
#include "fvModels.H"
#include "pisoControl.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "BoomsmaKeff.H"

#include <cmath>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "createFvModels.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"

    // Initialise phase-change and Boomsma conductivity fields at t=0.
    #include "updateCommonHotWallTemperature.H"
    #include "correctPCMProperties.H"
    #include "updateThermalBoundaryGradients.H"

    const scalar totalVolume = gSum(mesh.V());
    const scalar initialBulkTemperature =
        gSum(mesh.V()*Tp.primitiveField())/totalVolume;
    const scalar fixedBulkTemperature =
        TbulkReference > SMALL ? TbulkReference : initialBulkTemperature;

    const fileName dimensionlessDir
    (
        runTime.path()/"postProcessing"/"dimensionlessTime"
    );
    autoPtr<OFstream> dimensionlessTimeFilePtr;

    if (Pstream::master())
    {
        mkDir(dimensionlessDir);
        dimensionlessTimeFilePtr.reset
        (
            new OFstream(dimensionlessDir/"dimensionlessTime.dat")
        );

        dimensionlessTimeFilePtr()
            << "# Paper-comparison values use fixed liquid PCM properties."
            << nl
            << "# Ste=CpPl*|Tw-Tb|/Lf; alphaP=kPl/(rhoPl*CpPl)." << nl
            << "# Time TbFixed Tw CpSte alphaP meltFraction Ste Fo tau"
            << endl;
    }

    #include "writeDimensionlessTime.H"

    Info<< "Boomsma conductivities at t=0:" << nl
        << "    kPe  = " << gAverage(kPe)  << " W/(m K)" << nl
        << "    kSe  = " << gAverage(kSe)  << " W/(m K)" << nl
        << "    kEff = " << gAverage(kEff) << " W/(m K)" << nl
        << endl;

    const label nThermalCorrectors
    (
        piso.dict().lookupOrDefault<label>("nThermalCorrectors", 8)
    );

    Info<< "\nStarting LTNEPCMMetalFoamZhaoFoam time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "CourantNo.H"
        #include "setDeltaT.H"

        for
        (
            label thermalCorr = 0;
            thermalCorr < nThermalCorrectors;
            ++thermalCorr
        )
        {
            #include "updateCommonHotWallTemperature.H"
            #include "correctPCMProperties.H"
            #include "updateThermalBoundaryGradients.H"
            #include "TpEqn.H"

            #include "updateCommonHotWallTemperature.H"
            #include "correctPCMProperties.H"
            #include "updateThermalBoundaryGradients.H"
            #include "TsEqn.H"
        }

        #include "updateCommonHotWallTemperature.H"
        #include "correctPCMProperties.H"
        #include "updateThermalBoundaryGradients.H"

        if
        (
            gMin(Tp) < 100.0 || gMax(Tp) > 1000.0
         || gMin(Ts) < 100.0 || gMax(Ts) > 1000.0
        )
        {
            FatalErrorInFunction
                << "Unphysical temperature detected:" << nl
                << "    Tp min/max = " << gMin(Tp) << " " << gMax(Tp) << nl
                << "    Ts min/max = " << gMin(Ts) << " " << gMax(Ts) << nl
                << "Reduce maxDeltaT and inspect the thermal boundary "
                << "conditions and material properties."
                << exit(FatalError);
        }

        if (solveFlow)
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        runTime.write();

        if (runTime.writeTime())
        {
            #include "writeDimensionlessTime.H"
        }

        const volScalarField magU("magU", mag(U));

        Info<< "kPe=" << gAverage(kPe)
            << " kSe=" << gAverage(kSe)
            << " kEff=" << gAverage(kEff)
            << nl
            << "RePCM min/max=" << gMin(RePCM) << " " << gMax(RePCM)
            << "  PrPCM min/max=" << gMin(PrPCM) << " " << gMax(PrPCM)
            << nl
            << "|U| max=" << gMax(magU)
            << "  mushyCoeff min/max="
            << gMin(mushyCoeff) << " " << gMax(mushyCoeff)
            << nl
            << "hps min/max=" << gMin(hps) << " " << gMax(hps)
            << "  Tp min/max=" << gMin(Tp) << " " << gMax(Tp)
            << "  Ts min/max=" << gMin(Ts) << " " << gMax(Ts)
            << nl
            << "ExecutionTime = " << runTime.elapsedCpuTime()
            << " s  ClockTime = " << runTime.elapsedClockTime()
            << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
