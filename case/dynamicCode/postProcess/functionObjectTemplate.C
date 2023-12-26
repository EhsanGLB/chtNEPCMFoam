/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "functionObjectTemplate.H"
#include "functionObject.H"
#include "foamTime.H"
#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(postProcessFunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& postProcessFunctionObject::obr() const
{
    return obr_;
}


const fvMesh& postProcessFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

postProcessFunctionObject::postProcessFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

postProcessFunctionObject::~postProcessFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void postProcessFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read postProcess sha1: 3d452ce0e764058c4bf076968f631afcf4cad609\n";
    }

//{{{ begin code
    
//}}} end code
}


void postProcessFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute postProcess sha1: 3d452ce0e764058c4bf076968f631afcf4cad609\n";
    }

//{{{ begin code
    
//}}} end code
}


void postProcessFunctionObject::end()
{
    if (false)
    {
        Info<<"end postProcess sha1: 3d452ce0e764058c4bf076968f631afcf4cad609\n";
    }

//{{{ begin code
    
//}}} end code
}


void postProcessFunctionObject::timeSet()
{
    if (false)
    {
        Info<<"timeSet postProcess sha1: 3d452ce0e764058c4bf076968f631afcf4cad609\n";
    }

//{{{ begin codeTime
    
//}}} end code
}


void postProcessFunctionObject::write()
{
    if (false)
    {
        Info<<"write postProcess sha1: 3d452ce0e764058c4bf076968f631afcf4cad609\n";
    }

//{{{ begin code
    #line 78 "/home/ehsan/Desktop/chtNEPCMFoam/case/system/controlDict::functions::postProcess"
std::ofstream output;
            output.open("postProcess.txt",std::ofstream::app);

            //- dimensionless numbers
            scalar rhobf = 996.5;
            scalar mubf = 8.5e-4;
            scalar kappabf = 0.628;
            scalar Cpbf = 4181.0;
            scalar betabf = 2.0e-5;
            scalar nubf = mubf/rhobf;
            scalar alphabf = kappabf/(rhobf*Cpbf);
            scalar deltaTemp = 10.0;
            scalar charL = 1.0;
            scalar g = 6.43e-5;

            scalar Prbf = nubf/alphabf;
            scalar Grbf = (g*betabf*deltaTemp*pow(charL,3)) / (pow(nubf,2));
            scalar Rabf = (g*betabf*deltaTemp*pow(charL,3)) / (nubf*alphabf);


            //- lookup fields
            const volScalarField& T = mesh().lookupObject<volScalarField>("T");
            const volScalarField& Kappa = mesh().lookupObject<volScalarField>("Kappa");
            surfaceScalarField gradT = fvc::snGrad(T);

            //- find id and compute Nusselt number
            label leftID = T.mesh().boundaryMesh().findPatchID("left");
            scalar leftArea = sum(T.mesh().magSf().boundaryField()[leftID]);
            scalar leftHeatFluxAvg = sum(Kappa.boundaryField()[leftID]*gradT.boundaryField()[leftID]*T.mesh().magSf().boundaryField()[leftID])/leftArea;
            scalar lefthAvg = leftHeatFluxAvg/deltaTemp;
            scalar NusseltAvg = (lefthAvg*charL)/(kappabf);
            Info << "NusseltAvg = " << NusseltAvg << "\n";

            //- save data
            output << NusseltAvg << "\n";
            output.close();
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

