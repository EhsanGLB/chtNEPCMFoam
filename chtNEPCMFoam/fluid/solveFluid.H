
if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/UEqnFC.H"
    #include "forcedConvection/TEqnFC.H"

// --- PISO loop
    for (int corr = 0; corr < nCorr; corr++)
    {
        #include "forcedConvection/pEqnFC.H"
    }
}


if(typeFluid[i] == "naturalConvection")
{
    #include "naturalConvection/UEqnNC.H"
    #include "naturalConvection/TEqnNC.H"

// --- PISO loop
    for (int corr = 0; corr < nCorr; corr++)
    {
        #include "naturalConvection/pEqnNC.H"
    }
}


if(typeFluid[i] == "PCM")
{
    #include "PCM/UEqnPCM.H"
    #include "PCM/TEqnPCM.H"

// --- PISO loop
    for (int corr = 0; corr < nCorr; corr++)
    {
        #include "PCM/pEqnPCM.H"
    }
}


turb.correct();

