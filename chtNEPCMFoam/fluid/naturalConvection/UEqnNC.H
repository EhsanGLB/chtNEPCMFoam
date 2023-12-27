    Info << "    Type: Natural Convection" << endl;

    const dimensionedVector& g = gFluid[i];
    dimensionedScalar& TRef = TRefFluid[i];
    volScalarField& rhok = rhokFluid[i];


    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turb.divDevReff()
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
