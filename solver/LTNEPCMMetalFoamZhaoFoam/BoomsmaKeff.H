/*---------------------------------------------------------------------------*\
  Boomsma-Poulikakos tetrakaidecahedron effective conductivity model
  (paper Eqs. 13-19).

  ks      metal-foam skeleton conductivity [W/(m K)]
  kf      PCM conductivity [W/(m K)]
  epsilon metal-foam porosity [-]
  e       geometric constant = 0.339

  The function returns k_eff. The LTNE phase conductivities are formed as:
      kPe = k_eff(ks = 0,     kf = kPCM)
      kSe = k_eff(ks = kFoam, kf = 0)
\*---------------------------------------------------------------------------*/

#ifndef BoomsmaKeff_H
#define BoomsmaKeff_H

#include "scalar.H"
#include "mathematicalConstants.H"
#include <cmath>

inline Foam::scalar boomsmaKeff
(
    const Foam::scalar ks,
    const Foam::scalar kf,
    const Foam::scalar epsilon,
    const Foam::scalar e = 0.339
)
{
    using Foam::scalar;

    const scalar sqrt2 = std::sqrt(2.0);
    const scalar pi = Foam::constant::mathematical::pi;

    // Paper: d = sqrt{ sqrt(2)[2 - 5e^3 sqrt(2)/8 - 2 epsilon]
    //                    / [pi(3 - 4e sqrt(2) - e)] }
    const scalar dNumerator =
        sqrt2*(2.0 - 0.625*e*e*e*sqrt2 - 2.0*epsilon);
    const scalar dDenominator = pi*(3.0 - 4.0*e*sqrt2 - e);

    if
    (
        epsilon <= Foam::SMALL
     || epsilon >= 1.0 - Foam::SMALL
     || dNumerator <= Foam::SMALL
     || dDenominator <= Foam::SMALL
    )
    {
        FatalErrorInFunction
            << "Invalid Boomsma geometry: epsilon=" << epsilon
            << ", dNumerator=" << dNumerator
            << ", dDenominator=" << dDenominator
            << exit(Foam::FatalError);
    }

    const scalar d = std::sqrt(dNumerator/dDenominator);
    const scalar e2 = e*e;
    const scalar d2 = d*d;

    const scalar denA =
        (2.0*e2 + pi*d*(1.0 - e))*ks
      + (4.0 - 2.0*e2 - pi*d*(1.0 - e))*kf;

    const scalar denB =
        (e - 2.0*d)*e2*ks
      + (2.0*e - 4.0*d - (e - 2.0*d)*e2)*kf;

    const scalar denC =
        2.0*pi*d2*(1.0 - 2.0*e*sqrt2)*ks
      + 2.0
       *(
            sqrt2 - 2.0*e
          - pi*d2*(1.0 - e*sqrt2)
        )*kf;

    const scalar denD = e2*ks + (4.0 - e2)*kf;

    if
    (
        Foam::mag(denA) <= Foam::SMALL
     || Foam::mag(denB) <= Foam::SMALL
     || Foam::mag(denC) <= Foam::SMALL
     || Foam::mag(denD) <= Foam::SMALL
    )
    {
        FatalErrorInFunction
            << "Zero denominator in Boomsma conductivity model"
            << exit(Foam::FatalError);
    }

    const scalar RA = 4.0*d/denA;
    const scalar RB = Foam::sqr(e - 2.0*d)/denB;
    const scalar RC = Foam::sqr(sqrt2 - 2.0*e)/denC;
    const scalar RD = 2.0*e/denD;

    const scalar sumR = RA + RB + RC + RD;

    if (Foam::mag(sumR) <= Foam::SMALL)
    {
        FatalErrorInFunction
            << "Zero total resistance in Boomsma conductivity model"
            << exit(Foam::FatalError);
    }

    const scalar keff = sqrt2/(2.0*sumR);

    if (keff <= Foam::SMALL)
    {
        FatalErrorInFunction
            << "Non-positive Boomsma effective conductivity: " << keff
            << exit(Foam::FatalError);
    }

    return keff;
}

#endif

// ************************************************************************* //
