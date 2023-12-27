    Info << "    Type: PCM" << endl;

    const dimensionedVector& g = gFluid[i];
    dimensionedScalar& TRef = TRefFluid[i];
    volScalarField& rhok = rhokFluid[i];
    volScalarField& gamma = gammaFluid[i];
    dimensionedScalar& hm = hmFluid[i];
    dimensionedScalar& Tm = TmFluid[i];
    dimensionedScalar& Tmr = TmrFluid[i];
    dimensionedScalar& Cu = CuFluid[i];
    dimensionedScalar& b = bFluid[i];

    // Solve the Momentum equation
    volScalarField gammaV ("gammaV", Cu * pow( ( 1 - gamma ) , 2 ) / ( pow( gamma , 3 ) + b));
    Info << "max(p): " << max(p) << endl;


    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turb.divDevReff()
      + (1/rho)*fvm::SuSp(gammaV, U)
    );

    UEqn().relax();

    if (momentumPredictor)
    {
         solve
         (
            UEqn()
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rhok)*(g & mesh.Sf())
              - fvc::snGrad(p)*mesh.magSf()
            )
        );
    }
