    // Initialise fluid field pointer lists
    PtrList<volScalarField> TFluid(fluidRegions.size());
    PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> CpFluid(fluidRegions.size());
    PtrList<volScalarField> kappaFluid(fluidRegions.size());
    PtrList<volScalarField> muFluid(fluidRegions.size());
    PtrList<volScalarField> betaFluid(fluidRegions.size());
    PtrList<volScalarField> KappaFluid(fluidRegions.size());
    PtrList<word> typeFluid(fluidRegions.size());
    PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());
    PtrList<singlePhaseTransportModel> laminarTransport(fluidRegions.size());
    PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());





    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to TFluid\n" << endl;
        TFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to pFluid\n" << endl;
        pFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rho", dimDensity, SMALL)
            )
        );

        Info<< "    Adding to CpFluid\n" << endl;
        CpFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Cp",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL)
            )
        );

        Info<< "    Adding to kappaFluid\n" << endl;
        kappaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappa", dimThermalConductivity, SMALL)
            )
        );

        Info<< "    Adding to muFluid\n" << endl;
        muFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "mu",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("mu", dimMass/dimLength/dimTime, SMALL)
            )
        );

        Info<< "    Adding to betaFluid\n" << endl;
        betaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "beta",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("beta", dimless/dimTemperature, SMALL)
            )
        );

        Info<< "    Adding to KappaFluid\n" << endl;
        KappaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Kappa",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Kappa", dimThermalConductivity, SMALL)
            )
        );

        Info<< "    Adding to laminarTransport\n" << endl;
        laminarTransport.set
        (
            i,
            new singlePhaseTransportModel
            (
                UFluid[i],
                phiFluid[i]
            )
        );

        typeFluid.set
        (
            i,
            new word (laminarTransport[i].lookup("type"))
        );

        PrtFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Prt"))
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<incompressible::turbulenceModel>
            (
                incompressible::turbulenceModel::New
                (
                    UFluid[i],
                    phiFluid[i],
                    laminarTransport[i]
                )
            ).ptr()
        );
    }



// base fluid

    PtrList<dimensionedScalar> rhofFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappafFluid(fluidRegions.size());
    PtrList<dimensionedScalar> mufFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        rhofFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("rho"))
        );

        CpfFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Cp"))
        );

        kappafFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("kappa"))
        );

        mufFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("mu"))
        );
    }


// conditions

    // forced convection
    PtrList<volScalarField> rAUFluid(fluidRegions.size());
    PtrList<pimpleControl> pimpleFluid(fluidRegions.size());


    // natural convection and PCM
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betafFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TRefFluid(fluidRegions.size());
    PtrList<volScalarField> rhokFluid(fluidRegions.size());


    // PCM
    PtrList<volScalarField> gammaFluid(fluidRegions.size());
    PtrList<dimensionedScalar> hmFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmrFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CuFluid(fluidRegions.size());
    PtrList<dimensionedScalar> bFluid(fluidRegions.size());


    forAll(fluidRegions, i)
    {
      // forced convection
      if(typeFluid[i] == "forcedConvection")
      {
        Info<< "Reading field rAU if present\n" << endl;
        rAUFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rAU",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                runTime.deltaT(),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to pimpleFluid\n" << endl;
        pimpleFluid.set
        (
            i,
            new pimpleControl
            (
                fluidRegions[i],
                "PIMPLE"
            )
        );
      }


      // natural convection and PCM
      if(typeFluid[i] == "naturalConvection" || "PCM")
      {
        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        if(laminarTransport[i].found("beta"))
        {
            betafFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("beta"))
            );
        }

        if(laminarTransport[i].found("TRef"))
        {
            TRefFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("TRef"))
            );
        }

        Info<< "    Adding to rhokFluid\n" << endl;
        rhokFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhok",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhok", dimless, SMALL)
            )
        );
      }



      // PCM
      if(typeFluid[i] == "PCM")
      {
        Info<< "    Adding to gammaFluid\n" << endl;
        gammaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "gamma",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("gamma", dimless, SMALL)
            )
        );

        if(laminarTransport[i].found("hm"))
        {
            hmFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("hm"))
            );
        }

        if(laminarTransport[i].found("Tm"))
        {
            TmFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("Tm"))
            );
        }

        if(laminarTransport[i].found("Tmr"))
        {
            TmrFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("Tmr"))
            );
        }

        if(laminarTransport[i].found("Cu"))
        {
            CuFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("Cu"))
            );
        }

        if(laminarTransport[i].found("b"))
        {
            bFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("b"))
            );
        }
      }
    }

///////////NEPCM
    PtrList<IOdictionary> NEPCMProperties(solidRegions.size());
    PtrList<volScalarField> alphaFluid(fluidRegions.size());
    PtrList<volScalarField> rhoNEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> fNEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> NEPCM_moltenFluid(fluidRegions.size());
    PtrList<volScalarField> CpNEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> betaNEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> DBNEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> DTNEPCMFluid(fluidRegions.size());
    PtrList<word> alphaEqActivationFluid(fluidRegions.size());
    PtrList<dimensionedScalar> lFluid(fluidRegions.size());
    PtrList<dimensionedScalar> lambdaFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappaNEPCMFluid(fluidRegions.size());
    //-core
    PtrList<dimensionedScalar> rhocFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpcFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappacFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TMrFluid(fluidRegions.size());
    PtrList<dimensionedScalar> hcFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betacFluid(fluidRegions.size());
    PtrList<dimensionedScalar> dcFluid(fluidRegions.size());
    PtrList<volScalarField> CpPCMFluid(fluidRegions.size());
    //-shell
    PtrList<dimensionedScalar> rhosFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpsFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappasFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betasFluid(fluidRegions.size());
    PtrList<dimensionedScalar> dsFluid(fluidRegions.size());
   

    dimensionedScalar kB("kB", dimensionSet(1, 2, -2, -1, 0, 0, 0), 1.38064e-23);


    forAll(fluidRegions, i)
    {
        Info<< "    Adding to NEPCMProperties\n" << endl;
        NEPCMProperties.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "NEPCMProperties",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        //-NEPCM
        if (NEPCMProperties[i].headerOk())
        {
            Info<< "    Adding to alphaFluid\n" << endl;
            alphaFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "alpha",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i]
                )
            );

            Info<< "    Adding to rhoNEPCMFluid\n" << endl;
            rhoNEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "rhoNEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("rhoNEPCM", dimDensity, SMALL)
                )
            );

            Info<< "    Adding to fNEPCMFluid\n" << endl;
            fNEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "fNEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("fNEPCM", dimless, SMALL)
                )
            );

            Info<< "    Adding to NEPCM_moltenFluid\n" << endl;
            NEPCM_moltenFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "NEPCM_molten",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("NEPCM_molten", dimless, SMALL)
                )
            );

            Info<< "    Adding to CpNEPCMFluid\n" << endl;
            CpNEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "CpNEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("CpNEPCM", dimSpecificHeatCapacity, SMALL)
                )
            );

            Info<< "    Adding to betaNEPCMFluid\n" << endl;
            betaNEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betaNEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("betaNEPCM", dimless/dimTemperature, SMALL)
                )
            );

            alphaEqActivationFluid.set
            (
                i,
                new word (NEPCMProperties[i].lookup("alphaEqActivation"))
            );

            lFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("l"))
            );

            lambdaFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("lambda"))
            );

            //-core
            rhocFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("rhoc"))
            );

            CpcFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("Cpc"))
            );

            kappacFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("kappac"))
            );

            TfFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("Tf"))
            );

            TMrFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("TMr"))
            );

            hcFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("hc"))
            );

            if(NEPCMProperties[i].found("betac"))
            {
                betacFluid.set
                (
                    i,
                    new dimensionedScalar (NEPCMProperties[i].lookup("betac"))
                );
            }

            Info<< "    Adding to CpPCMFluid\n" << endl;
            CpPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "CpPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("CpPCM", dimSpecificHeatCapacity, SMALL)
                )
            );

            //-shell
            rhosFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("rhos"))
            );

            CpsFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("Cps"))
            );

            kappasFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("kappas"))
            );

            if(NEPCMProperties[i].found("betas"))
            {
                betasFluid.set
                (
                    i,
                    new dimensionedScalar (NEPCMProperties[i].lookup("betas"))
                );
            }

            dsFluid.set
            (
                i,
                new dimensionedScalar (NEPCMProperties[i].lookup("ds"))
            );

            kappaNEPCMFluid.set
            (
                i,
                new dimensionedScalar ("kappaNEPCM", dimThermalConductivity, SMALL)
            );

            dcFluid.set
            (
                i,
                new dimensionedScalar ("dc", dimLength, SMALL)
            );
        }

    }






