/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "exampleClassModel.H"
#include "vector.H"
#include "scalarList.H"
#include "labelList.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * /.
// This declares an external function call to the fortran function in LAPACK
extern "C" void dgesv_(const long *Np, const long *NRHSp, double *A, long *LDAp, double *IPIV, double *B, const long *LDBp, long *INFOp);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Declare a C++ function to interface with the LAPACK external function
static long dgesv(long N, long NRHS, double *A, long LDA, double *IPIV, double *B, long LDB)
{
    long info;
    dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &info);
    return info;
}

int exampleClassModel::readModelData()
{
    IFstream inputStream("system/exampleDict");

    if (inputStream.good())
    {
        dictionary dict(inputStream);
        Info << "Reading exampleDict" << endl;

        dict.lookup("N") >> nModes_;
    }
    else
    {
        FatalErrorIn
        (
            "exampleClassModel::readModeData()"
        )   << "Cannot open exampleDict"
            << exit(FatalError);
    }
    return 0;
}

tmp< Field<scalarField> > exampleClassModel::mapFluid2Example(const vectorField& fluidPts, const vectorField& faceCenters, vectorField n)
{

    tmp< Field<scalarField> > matching(new Field<scalarField>(2));
    scalarField& faceMatching = matching()[0];
    faceMatching.setSize(faceCenters.size(), -1);
    scalarField& pointMatching = matching()[1];
    pointMatching.setSize(fluidPts.size(), -1);

    fluidFaces2ExpMap_.setSize(faceCenters.size());
    fluidFaces2ExpWeights_.setSize(faceCenters.size());

    return matching;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct
exampleClassModel::exampleClassModel
(
    word name,
    bool uncoupled,
    scalar& dt,
    fileName timePath
)
:
    uncoupled_(uncoupled),
    timePath_(timePath),
    name_(name)
{
    readModelData();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::exampleClassModel::~exampleClassModel(){}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp< Field<scalarField> > exampleClassModel::map2fluidDomain
(
    const vectorField& fluidPts,
    const vectorField& fluidFaceCenters,
    vectorField& fluidFaceNormals
)
{
    return mapFluid2Example(fluidPts, fluidFaceCenters, fluidFaceNormals);
}

void exampleClassModel::updateState(scalarField& p, vector Ag, vector alpha)
{
}

void exampleClassModel::updateState(scalarField& p, vector Ag, vector alpha, scalar dt)
{
}

void exampleClassModel::updateState(scalarField& p, vector Ag, vector alpha, word timeName)
{
}

void exampleClassModel::updateState(scalarField& p, vector Ag, vector alpha, scalar dt, word timeName)
{
}

void exampleClassModel::advanceState(scalarField& p, vector Ag, vector alpha)
{
}

void exampleClassModel::advanceState(scalarField& p, vector Ag, vector alpha, scalar dt)
{
}


// Access

// Return const access to the current velocities at the fluid faces
tmp<vectorField> exampleClassModel::currentFaceVelocity()
{
    tmp<vectorField> tfaceVelocity(new vectorField(fluidFaces2ExpMap_.size(), vector::zero));
    vectorField& faceVelocity = tfaceVelocity();
    forAll(faceVelocity, ii)
    {
        for (int nn = 0; nn < nModes_; nn++)
        {
            vector w(vector::zero);
            forAll(fluidFaces2ExpMap_[ii], jj)
            {
                if (fluidFaces2ExpMap_[ii][jj] != -1)
                    for (int kk = 0; kk < 6; kk++)
                        w += fluidFaces2ExpWeights_[ii][6*jj+kk];
            }
            faceVelocity[ii] += w;
        }
    }

    return tfaceVelocity;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
