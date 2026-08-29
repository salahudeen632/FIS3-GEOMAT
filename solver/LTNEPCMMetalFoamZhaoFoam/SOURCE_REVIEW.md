# Reviewed source decisions

## Correct LTNE conduction and exchange terms

```cpp
// PCM
- fvm::laplacian(kPe, Tp)
+ fvm::Sp(hps*aPS, Tp)
== hps*aPS*Ts;

// foam
- fvm::laplacian(kSe, Ts)
+ fvm::Sp(hps*aPS, Ts)
== hps*aPS*Tp;
```

These terms conserve interfacial exchange: the PCM gain is the foam loss.

## No-flow limitation

The supplied paper gives `hps` as a Reynolds-number correlation but does not
provide a stationary (`Re=0`) closure. A pure no-flow use of Eq. (12) gives
`hps=0`. This solver exposes that choice explicitly through `hpsModel`.

## Reference velocity for Reynolds number

The code uses:

```text
RePCM = rhoP*mag(U)*dp/muP
```

where `U` is the volume-averaged/superficial velocity solved by the model.
The supplied document names `Re` but does not explicitly distinguish
superficial from interstitial velocity; keep this definition unless the source
for Eq. (12) specifies another convention.

## PCM apparent-heat-capacity transient

```cpp
fvm::Sp(epsilon*rhoCp/deltaT, Tp)
==
epsilon*rhoCp/deltaT*Tp.oldTime() + ...
```

This directly discretizes `epsilon*rho*Cp*dTp/dt` using the current Picard
iterate of apparent heat capacity. It avoids differentiating the sharply
temperature-dependent latent-heat peak, which can produce severe temperature
undershoots.

## Deliberately retained regularizer

Per the requested model choice, the mushy function remains:

```text
(1-b)^2/(b^3 + q), q = 1e-3
```
