Info<< "Reading g\n" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading pcmMetalFoamProperties\n" << endl;

IOdictionary pcmMetalFoamProperties
(
    IOobject
    (
        "pcmMetalFoamProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const Switch solveFlow
(
    pcmMetalFoamProperties.lookupOrDefault<Switch>("solveFlow", true)
);

const scalar epsilon
(
    readScalar(pcmMetalFoamProperties.lookup("epsilon"))
);

if (epsilon <= 0.0 || epsilon >= 1.0)
{
    FatalErrorInFunction
        << "epsilon must be strictly between 0 and 1. Supplied epsilon="
        << epsilon << exit(FatalError);
}

const dimensionedScalar rhoPs
(
    "rhoPs", dimDensity, pcmMetalFoamProperties.lookup("rhoPs")
);
const dimensionedScalar rhoPl
(
    "rhoPl", dimDensity, pcmMetalFoamProperties.lookup("rhoPl")
);
const dimensionedScalar CpPs
(
    "CpPs", dimSpecificHeatCapacity, pcmMetalFoamProperties.lookup("CpPs")
);
const dimensionedScalar CpPl
(
    "CpPl", dimSpecificHeatCapacity, pcmMetalFoamProperties.lookup("CpPl")
);
const dimensionedScalar kPs
(
    "kPs", dimPower/(dimLength*dimTemperature), pcmMetalFoamProperties.lookup("kPs")
);
const dimensionedScalar kPl
(
    "kPl", dimPower/(dimLength*dimTemperature), pcmMetalFoamProperties.lookup("kPl")
);
const dimensionedScalar muPl
(
    "muPl", dimDynamicViscosity, pcmMetalFoamProperties.lookup("muPl")
);
const dimensionedScalar Lf
(
    "Lf", dimEnergy/dimMass, pcmMetalFoamProperties.lookup("Lf")
);
const dimensionedScalar Tm
(
    "Tm", dimTemperature, pcmMetalFoamProperties.lookup("Tm")
);
const word phaseChangeDeltaTEntry
(
    pcmMetalFoamProperties.found("phaseChangeDeltaT")
  ? word("phaseChangeDeltaT")
  : word("deltaT")
);
const dimensionedScalar phaseChangeDeltaT
(
    "phaseChangeDeltaT",
    dimTemperature,
    pcmMetalFoamProperties.lookup(phaseChangeDeltaTEntry)
);
const dimensionedScalar beta
(
    "beta", dimless/dimTemperature, pcmMetalFoamProperties.lookup("beta")
);
const dimensionedScalar Tref
(
    "Tref", dimTemperature, pcmMetalFoamProperties.lookup("Tref")
);
const dimensionedScalar Tinf
(
    "Tinf", dimTemperature, pcmMetalFoamProperties.lookup("Tinf")
);
const dimensionedScalar Tw
(
    "Tw",
    dimTemperature,
    pcmMetalFoamProperties.lookupOrDefault<scalar>("Tw", Tinf.value())
);
const word hotWallThermalMode
(
    pcmMetalFoamProperties.lookupOrDefault<word>
    (
        "hotWallThermalMode",
        "none"
    )
);
const dimensionedScalar qHotWall
(
    "qHotWall",
    dimPower/dimArea,
    pcmMetalFoamProperties.lookupOrDefault<scalar>("qHotWall", 0.0)
);
const dimensionedScalar qBottom
(
    "qBottom", dimPower/dimArea, pcmMetalFoamProperties.lookup("qBottom")
);
const dimensionedScalar hLoss
(
    "hLoss",
    dimPower/(dimArea*dimTemperature),
    pcmMetalFoamProperties.lookup("hLoss")
);

const word hotWallPatchName
(
    pcmMetalFoamProperties.lookupOrDefault<word>("hotWallPatch", "leftWall")
);

if
(
    hotWallThermalMode != "fixedTemperature"
 && hotWallThermalMode != "heatFlux"
 && hotWallThermalMode != "none"
)
{
    FatalErrorInFunction
        << "hotWallThermalMode must be fixedTemperature, heatFlux, or none. "
        << "Supplied value: " << hotWallThermalMode
        << exit(FatalError);
}

if (hotWallThermalMode == "heatFlux" && qHotWall.value() <= SMALL)
{
    FatalErrorInFunction
        << "hotWallThermalMode is heatFlux, but qHotWall is not positive: "
        << qHotWall << exit(FatalError);
}

const scalar characteristicLength
(
    pcmMetalFoamProperties.lookupOrDefault<scalar>
    (
        "characteristicLength",
        max
        (
            mesh.bounds().span().x(),
            max(mesh.bounds().span().y(), mesh.bounds().span().z())
        )
    )
);

const scalar TbulkReference
(
    pcmMetalFoamProperties.lookupOrDefault<scalar>
    (
        "TbulkReference",
        0.0
    )
);

const dimensionedScalar rhoFoam
(
    "rhoFoam", dimDensity, pcmMetalFoamProperties.lookup("rhoFoam")
);
const dimensionedScalar CpFoam
(
    "CpFoam", dimSpecificHeatCapacity, pcmMetalFoamProperties.lookup("CpFoam")
);
const dimensionedScalar kFoam
(
    "kFoam", dimPower/(dimLength*dimTemperature), pcmMetalFoamProperties.lookup("kFoam")
);

const scalar PPI
(
    readScalar(pcmMetalFoamProperties.lookup("PPI"))
);
const scalar viscosityA
(
    readScalar(pcmMetalFoamProperties.lookup("viscosityA"))
);
const dimensionedScalar mushyCu
(
    "mushyCu", dimless/dimTime, pcmMetalFoamProperties.lookup("mushyCu")
);
const scalar q
(
    readScalar(pcmMetalFoamProperties.lookup("q"))
);
const Switch useViscosityDamping
(
    pcmMetalFoamProperties.lookupOrDefault<Switch>
    (
        "useViscosityDamping", true
    )
);
const word hpsModel
(
    pcmMetalFoamProperties.lookupOrDefault<word>
    (
        "hpsModel",
        (solveFlow ? word("Zukauskas") : word("constant"))
    )
);

// hpsConstant and hpsMin are entered as numerical values in W/(m2 K).
const dimensionedScalar hpsConstant
(
    "hpsConstant",
    dimPower/(dimArea*dimTemperature),
    pcmMetalFoamProperties.lookupOrDefault<scalar>("hpsConstant", 0.0)
);

const dimensionedScalar hpsMin
(
    "hpsMin",
    dimPower/(dimArea*dimTemperature),
    pcmMetalFoamProperties.lookupOrDefault<scalar>("hpsMin", 0.0)
);

if (phaseChangeDeltaT.value() <= SMALL)
{
    FatalErrorInFunction
        << "phaseChangeDeltaT must be > 0. Supplied phaseChangeDeltaT="
        << phaseChangeDeltaT
        << exit(FatalError);
}

if (q <= SMALL)
{
    FatalErrorInFunction
        << "q must be > 0. Supplied q=" << q
        << exit(FatalError);
}

if
(
    viscosityA < 0.0
 || mushyCu.value() < 0.0
 || hpsConstant.value() < 0.0
 || hpsMin.value() < 0.0
 || qBottom.value() < 0.0
 || hLoss.value() < 0.0
)
{
    FatalErrorInFunction
        << "viscosityA, mushyCu, hpsConstant, hpsMin, qBottom and hLoss "
        << "must be non-negative."
        << exit(FatalError);
}

if (useViscosityDamping && viscosityA > SMALL && mushyCu.value() > SMALL)
{
    WarningInFunction
        << "Both viscosity damping and the mushy momentum source are active. "
        << "This applies solidification resistance twice. For the standard "
        << "enthalpy-porosity source formulation, set useViscosityDamping "
        << "false and viscosityA 0." << nl;
}

if (hpsModel != "Zukauskas" && hpsModel != "constant")
{
    FatalErrorInFunction
        << "hpsModel must be either Zukauskas or constant. Supplied hpsModel="
        << hpsModel << exit(FatalError);
}

if (!solveFlow && hpsModel == "Zukauskas" && hpsMin.value() <= SMALL)
{
    WarningInFunction
        << "solveFlow is false and hpsModel is Zukauskas with hpsMin = 0. "
        << "This gives Re = 0 and therefore hps = 0, decoupling Tp and Ts. "
        << "Use hpsModel constant or prescribe a positive hpsMin for a "
        << "no-flow LTNE calculation." << nl;
}

if (PPI <= SMALL)
{
    FatalErrorInFunction
        << "PPI must be > 0. Supplied PPI=" << PPI
        << exit(FatalError);
}

if (characteristicLength <= SMALL)
{
    FatalErrorInFunction
        << "characteristicLength must be > 0. Supplied value="
        << characteristicLength << exit(FatalError);
}

// Metal-foam microstructure correlations from Eqs. (20)-(23).
const scalar dp = 0.0254/PPI;
const scalar dfOverDp =
    1.18*std::sqrt((1.0 - epsilon)/constant::mathematical::pi)
   /(1.0 - std::exp(-(1.0 - epsilon)/0.04));
const scalar df = dfOverDp*dp;

const dimensionedScalar K
(
    "K",
    dimArea,
    0.00073*sqr(dp)*std::pow(1.0 - epsilon, -0.224)
   *std::pow(dfOverDp, -1.11)
);

const dimensionedScalar CForch
(
    "CForch",
    dimless,
    // Reference corrected solver form.
    0.00212*sqr(dp)*std::pow(1.0 - epsilon, -0.132)
   *std::pow(dfOverDp, -1.63)
);

const dimensionedScalar aPS
(
    "aPS",
    dimless/dimLength,
    3.0*constant::mathematical::pi*df
   *(1.0 - std::exp(-(1.0 - epsilon)/0.004))
   /sqr(0.59*dp)
);

Info<< "Metal-foam correlations:" << nl
    << "    pore diameter dp       = " << dp << " m" << nl
    << "    ligament ratio df/dp   = " << dfOverDp << nl
    << "    permeability K         = " << K << nl
    << "    inertial coefficient C = " << CForch << nl
    << "    specific area aPS      = " << aPS << nl
    << "    hpsModel               = " << hpsModel << nl
    << "    hpsConstant            = " << hpsConstant << nl
    << "    hpsMin                 = " << hpsMin << nl
    << "    qBottom                = " << qBottom << nl
    << "    hLoss                  = " << hLoss << nl
    << "    Tinf                   = " << Tinf << nl
    << "    Tw                     = " << Tw << nl
    << "    phaseChangeDeltaT      = " << phaseChangeDeltaT << nl
    << "    hotWallThermalMode     = " << hotWallThermalMode << nl
    << "    qHotWall               = " << qHotWall << nl
    << "    solveFlow              = " << solveFlow << nl
    << endl;

Info<< "Creating gh\n" << endl;
volScalarField gh
(
    IOobject
    (
        "gh",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    g & mesh.C()
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Tp (PCM temperature)\n" << endl;
volScalarField Tp
(
    IOobject
    (
        "Tp",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ts (metal-foam temperature)\n" << endl;
volScalarField Ts
(
    IOobject
    (
        "Ts",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

volScalarField alpha
(
    IOobject
    (
        "alphaPCM",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField deltaDirac
(
    IOobject
    (
        "deltaDirac",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTemperature, 0.0)
);

volScalarField rhoP
(
    IOobject("rhoP", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    rhoPs
);

volScalarField p
(
    IOobject("p", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    p_rgh + rhoP*gh
);

volScalarField CpP
(
    IOobject("CpP", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    CpPs
);

volScalarField CpSensible
(
    IOobject
    (
        "CpSensible",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    CpPs
);

volScalarField rhoCp
(
    IOobject("rhoCp", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    rhoPs*CpPs
);

volScalarField kP
(
    IOobject("kP", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    kPs
);

volScalarField muP
(
    IOobject("muP", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    muPl
);

volScalarField kPe
(
    IOobject("kPe", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    kPs
);

volScalarField kSe
(
    IOobject("kSe", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    kFoam
);

volScalarField kEff
(
    IOobject("kEff", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    kPs
);

volScalarField hps
(
    IOobject("hps", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    hpsMin
);

volScalarField RePCM
(
    IOobject("RePCM", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField PrPCM
(
    IOobject("PrPCM", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField mushyCoeff
(
    IOobject("mushyCoeff", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p_rgh, piso.dict(), pRefCell, pRefValue);
mesh.schemes().setFluxRequired(p_rgh.name());

const label hotWallPatchi =
    hotWallThermalMode == "none"
  ? -1
  : mesh.boundaryMesh().findPatchID(hotWallPatchName);

if (hotWallThermalMode != "none" && hotWallPatchi < 0)
{
    FatalErrorInFunction
        << "Cannot find hot-wall patch " << hotWallPatchName
        << exit(FatalError);
}
