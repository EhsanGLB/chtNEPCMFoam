    const fvMesh& mesh = fluidRegions[i];

    volScalarField& T = TFluid[i];
    volScalarField& p = pFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& rho = rhoFluid[i];
    volScalarField& Cp = CpFluid[i];
    volScalarField& kappa = kappaFluid[i];
    volScalarField& mu = muFluid[i];
    volScalarField& beta = betaFluid[i];
    volScalarField& Kappa = KappaFluid[i];
    word& type = typeFluid[i];
    dimensionedScalar& Prt = PrtFluid[i];
    incompressible::turbulenceModel& turb = turbulence[i];

    mesh.schemesDict().setFluxRequired(p.name());

   //-base fluid
    dimensionedScalar& rhof = rhofFluid[i];
    dimensionedScalar& Cpf = CpfFluid[i];
    dimensionedScalar& kappaf = kappafFluid[i];
    dimensionedScalar& muf = mufFluid[i];

    if (NEPCMProperties[i].headerOk())
    {
        //-NEPCM
        volScalarField& alpha = alphaFluid[i];
        volScalarField& rhoNEPCM = rhoNEPCMFluid[i];
        volScalarField& fNEPCM = fNEPCMFluid[i];
        volScalarField& NEPCM_molten = NEPCM_moltenFluid[i];
        volScalarField& CpNEPCM = CpNEPCMFluid[i];
        word& alphaEqActivation = alphaEqActivationFluid[i];
        dimensionedScalar& l = lFluid[i];
        dimensionedScalar& lambda = lambdaFluid[i];
        dimensionedScalar& kappaNEPCM = kappaNEPCMFluid[i];
        //-core
        dimensionedScalar& rhoc = rhocFluid[i];
        dimensionedScalar& Cpc = CpcFluid[i];
        dimensionedScalar& kappac = kappacFluid[i];
        dimensionedScalar& Tf = TfFluid[i];
        dimensionedScalar& TMr = TMrFluid[i];
        dimensionedScalar& hc = hcFluid[i];
        dimensionedScalar& dc = dcFluid[i];
        volScalarField& CpPCM = CpPCMFluid[i];
        //-shell
        dimensionedScalar& rhos = rhosFluid[i];
        dimensionedScalar& Cps = CpsFluid[i];
        dimensionedScalar& kappas = kappasFluid[i];
        dimensionedScalar& ds = dsFluid[i];
        //-NEPCM
        dc = ds*pow( rhos/(rhos+l*rhoc), scalar(0.3333) );
        kappaNEPCM = ds/( (dc/kappac) + ((ds-dc)/kappas) );

        forAll(fNEPCM, j)
        {
            if (T[j] < scalar(Tf.value()-TMr.value()/2) )
            {
                fNEPCM[j] = scalar(0.0);
                NEPCM_molten[j] = scalar(0.0);
            }
            else if ( scalar(Tf.value()-TMr.value()/2) <= T[j] && T[j] <= scalar(Tf.value()+TMr.value()/2) )
            {
                fNEPCM[j] = scalar(1.0);
                NEPCM_molten[j] = scalar(0.5);
            }
            else if (T[j] > scalar(Tf.value()+TMr.value()/2) )
            {
                fNEPCM[j] = scalar(0.0);
                NEPCM_molten[j] = scalar(1.0);
            }
        }

        rhoNEPCM = ((1+l)*rhoc*rhos)/(rhos+l*rhoc);
        CpPCM = Cpc+( (3.14/2)*(hc/TMr-Cpc)*sin( (3.14/TMr)*(T-(Tf-TMr/2)) ) )*fNEPCM;
        CpNEPCM = ((CpPCM+l*Cps)*rhoc*rhos)/((rhos+l*rhoc)*rhoNEPCM);

        //-mixture
        rho = alpha*rhoNEPCM+(1-alpha)*rhof;
        Cp = (alpha*rhoNEPCM*CpNEPCM+(1-alpha)*rhof*Cpf)/rho;
        kappa = ((kappaNEPCM+(lambda-1)*kappaf-(lambda-1)*(kappaf-kappaNEPCM)*alpha)/(kappaNEPCM+(lambda-1)*kappaf+(kappaf-kappaNEPCM)*alpha))*kappaf;
        mu = (muf/pow((1.0-alpha),2.5));

        if(type == "naturalConvection" || "PCM")
        {
            dimensionedScalar& betaf = betafFluid[i];
            volScalarField& betaNEPCM = betaNEPCMFluid[i];
            dimensionedScalar& betac = betacFluid[i];
            dimensionedScalar& betas = betasFluid[i];
            betaNEPCM = betac + ( (betas-betac)/2 )*(1 - l*rhos/rhoc);
            beta = (alpha*rhoNEPCM*betaNEPCM+(1-alpha)*rhof*betaf)/rho;
        }
    }
    else
    {
        rho = rhof;
        Cp = Cpf;
        kappa = kappaf;
        mu = muf;

        if(type == "naturalConvection" || "PCM")
        {
            dimensionedScalar& betaf = betafFluid[i];
            beta = betaf;
        }
    }
