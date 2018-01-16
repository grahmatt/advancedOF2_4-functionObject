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

#include "advanceFunctionObject.H"
#include "volFields.H"
#include "pointFields.H"
#include "pointPatchField.H"
#include "uniformDimensionedFields.H"
#include "dictionary.H"
#include "Time.H"
#include "SHA1Digest.H"
#include "stringOps.H"
#include "ListListOps.H"
#include "addToRunTimeSelectionTable.H"
#include <iostream>
#include "forces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(advanceFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        advanceFunctionObject,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::advanceFunctionObject::description() const
{
    return "functionObject " + name();
}


const Foam::dictionary& Foam::advanceFunctionObject::exampleDict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::advanceFunctionObject::advanceFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict,
    bool readNow
)
:
    functionObject(name),
    time_(time),
    dict_(dict),
    exampleDict_
    (
        IOobject
        (
            "exampleDict",
            time.caseSystem(),
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    obr_
    (
        time.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault<word>("region",polyMesh::defaultRegion)
        )
    ),
    patchIDs_(0,-1),
    adjustTimeStep_(time.controlDict().lookupOrDefault<Switch>("adjustTimeStep",false)),
    pName_(exampleDict_.lookupOrDefault<word>("p","p")),
    exampleModel_
    (
        exampleDict_.lookupOrDefault<word>("modelName",""),
        dict.lookupOrDefault<Switch>("uncoupledModes",Switch::TRUE),
        time.deltaT().value(),
        (Pstream::parRun() ? (time.path()/"..") : time.path())
    ),
    faceStarts_(1,0),
    pointStarts_(1,0),
    patchFaceStarts_(1,0),
    patchPointStarts_(1,0),
    faceCenters_(0,vector::zero),
    fluidPts_(0,vector::zero),
    outputInterval_(dict.lookupOrDefault<label>("outputInterval",1)),
    curTimeIndex_(-1)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::advanceFunctionObject::~advanceFunctionObject() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::advanceFunctionObject::start()
{
    if (!obr_.foundObject<volScalarField>(pName_))
        return false;

    const volScalarField& pressure = obr_.lookupObject<volScalarField>(pName_);

    // obtain fvMesh to obtain fluid faces
    const fvMesh& mesh = pressure.mesh();

    int nPoints = 0, nFaces = 0;
    if (patchIDs_.size() == 0)
    {
        wordList patchNames(exampleDict_.lookup("patches"));
        patchIDs_.setSize(patchNames.size(),-1);
        patchFaceStarts_.setSize(patchNames.size()+1,0);
        patchPointStarts_.setSize(patchNames.size()+1,0);
        forAll(patchIDs_, ii)
        {
            patchIDs_[ii] = mesh.boundary().findPatchID(patchNames[ii]);
            nFaces += mesh.boundary()[patchIDs_[ii]].Cf().size();
            nPoints += mesh.boundaryMesh()[patchIDs_[ii]].localPoints().size();
            patchFaceStarts_[ii+1] = nFaces;
            patchPointStarts_[ii+1] = nPoints;
        }
    }

    List<vectorField> faceCentersList(Pstream::nProcs());
    List<vectorField> fluidPtsList(Pstream::nProcs());
    List<vectorField> nList(Pstream::nProcs());
    vectorField combPatchFaceCenters(nFaces,vector::zero);
    vectorField combPatchFluidPts(nPoints,vector::zero);
    vectorField combPatchNormals(nFaces,vector::zero);
    forAll(patchIDs_, ii)
    {
        for (int jj = patchFaceStarts_[ii]; jj < patchFaceStarts_[ii+1]; jj++)
        {
            combPatchFaceCenters[jj] = mesh.boundary()[patchIDs_[ii]].Cf()[jj-patchFaceStarts_[ii]];
            combPatchNormals[jj] = mesh.boundary()[patchIDs_[ii]].nf()()[jj-patchFaceStarts_[ii]];
        }
        for (int jj = patchPointStarts_[ii]; jj < patchPointStarts_[ii+1]; jj++)
            combPatchFluidPts[jj] = mesh.boundaryMesh()[patchIDs_[ii]].localPoints()[jj-patchPointStarts_[ii]];
    }
    faceCentersList[Pstream::myProcNo()] = combPatchFaceCenters;
    fluidPtsList[Pstream::myProcNo()] = combPatchFluidPts;
    nList[Pstream::myProcNo()] = combPatchNormals;
    Pstream::gatherList(faceCentersList);
    Pstream::gatherList(fluidPtsList);
    Pstream::gatherList(nList);

    PstreamBuffers buf(Pstream::defaultCommsType);
    if (Pstream::master())
    {
        // find start of each processor in agglomerated lists
        faceStarts_.setSize(Pstream::nProcs()+1,0);
        pointStarts_.setSize(Pstream::nProcs()+1,0);
        for (int ii = 0; ii < Pstream::nProcs(); ii++)
        {
            faceStarts_[ii+1] = faceStarts_[ii]+faceCentersList[ii].size();
            pointStarts_[ii+1] = pointStarts_[ii]+fluidPtsList[ii].size();
        }

        // agglomerate lists from each processor into one
        faceCenters_ =
            ListListOps::combine<vectorField >
            (
                faceCentersList,
                accessOp<vectorField >()
            );
        fluidPts_ =
            ListListOps::combine<vectorField >
            (
                fluidPtsList,
                accessOp<vectorField >()
            );
        vectorField n
        (
            ListListOps::combine<vectorField >
            (
                nList,
                accessOp<vectorField >()
            )
        );

        for (int ii = 0; ii < Pstream::nProcs(); ii++)
        {
            if (ii != Pstream::myProcNo())
            {
                vectorField procOut(faceCenters_.size() + fluidPts_.size() + n.size() + 1);
                vector sizes(faceCenters_.size(), fluidPts_.size(), n.size());
                procOut[0] = sizes;
                for (int jj = 0; jj < faceCenters_.size(); jj++)
                    procOut[jj+1] = faceCenters_[jj];
                for (int jj = 0; jj < fluidPts_.size(); jj++)
                    procOut[jj+faceCenters_.size()+1] = fluidPts_[jj];
                for (int jj = 0; jj < n.size(); jj++)
                    procOut[jj+faceCenters_.size()+fluidPts_.size()+1] = n[jj];

                UOPstream toNeighbor(ii, buf);
                toNeighbor << procOut;
            }
        }

        // send fluid grid information example class model
        exampleModel_.map2fluidDomain(fluidPts_, faceCenters_, n);

        buf.finishedSends();
    }
    else
    {
        vectorField procIn;
        buf.finishedSends();
        UIPstream fromMaster(Pstream::masterNo(), buf);
        fromMaster >> procIn;

        faceCenters_.resize(static_cast<int>(procIn[0].x()+0.1));
        fluidPts_.resize(static_cast<int>(procIn[0].y()+0.1));
        vectorField n(static_cast<int>(procIn[0].z()+0.1));
        for (int jj = 0; jj < faceCenters_.size(); jj++)
            faceCenters_[jj] = procIn[jj+1];
        for (int jj = 0; jj < fluidPts_.size(); jj++)
            fluidPts_[jj] = procIn[jj+faceCenters_.size()+1];
        for (int jj = 0; jj < n.size(); jj++)
            n[jj] = procIn[jj+faceCenters_.size()+fluidPts_.size()+1];
    }

    return true;
}


bool Foam::advanceFunctionObject::execute(const bool forceWrite)
{
    if (!obr_.foundObject<volScalarField>(pName_))
        return false;

    List<scalarField> pList(Pstream::nProcs());

    if (patchIDs_.size() > 0)
    {
        const volScalarField& pressure = obr_.lookupObject<volScalarField>(pName_);

        scalarField p(patchFaceStarts_[patchIDs_.size()],0.0);
        forAll(patchIDs_, ii)
        {
            scalarField patchp = pressure.boundaryField()[patchIDs_[ii]];
            for (int jj = patchFaceStarts_[ii]; jj < patchFaceStarts_[ii+1]; jj++)
                p[jj] = patchp[jj-patchFaceStarts_[ii]];
        }

        pList[Pstream::myProcNo()] = p;
        Pstream::gatherList(pList);
    }
    else
    {
        scalarField tmp(0);
        pList[Pstream::myProcNo()] = tmp;
        Pstream::gatherList(pList);
    }

    PstreamBuffers buf(Pstream::defaultCommsType);
    scalarField p;
    if (Pstream::master())
    {
        // Combine pressure values into single field
        p = ListListOps::combine<scalarField >
            (
                pList,
                accessOp<scalarField >()
            );
        if (p.size() != faceStarts_[Pstream::nProcs()])
        {
            FatalErrorIn
            (
                "advanceFunctionObject::execute(const bool)"
            )   << "combined pressure and faceCenters are not same size (" << p.size() << " versus " << faceStarts_[Pstream::nProcs()] << ") on processor " << Pstream::myProcNo()
                << exit(FatalError);
        }

        for (int ii = 0; ii < Pstream::nProcs(); ii++)
        {
            if (ii != Pstream::myProcNo())
            {
                UOPstream toNeighbor(ii, buf);
                toNeighbor << p;
            }
        }
        buf.finishedSends();
    }
    else
    {
        buf.finishedSends();
        UIPstream fromMaster(Pstream::masterNo(), buf);
        fromMaster >> p;
    }

    // initialize scattered values
    List<vectorField> outputVectorList(Pstream::nProcs());
    if (Pstream::myProcNo() < nExampleProc_)
    {
        vector g(vector::zero);
        if (obr_.foundObject<uniformDimensionedVectorField>("g"))
            g = obr_.lookupObject<uniformDimensionedVectorField>("g").value();

        vector a(vector::zero), alpha(vector::zero);

        if (curTimeIndex_ != time_.timeIndex()) // advance state
        {
            if (adjustTimeStep_)
                exampleModel_.advanceState(p, a-g, alpha, time_.deltaT().value());
            else
                exampleModel_.advanceState(p, a-g, alpha);

            curTimeIndex_ = time_.timeIndex();
        }
        else // update state for iteration
        {

            if (adjustTimeStep_)
            {
                if (curTimeIndex_%outputInterval_ == 0)
                    exampleModel_.updateState(p, a-g, alpha, time_.deltaT().value(), time_.timeName());
                else
                    exampleModel_.updateState(p, a-g, alpha, time_.deltaT().value());
            }
            else
            {
                if (curTimeIndex_%outputInterval_ == 0)
                    exampleModel_.updateState(p, a-g, alpha, time_.timeName());
                else
                    exampleModel_.updateState(p, a-g, alpha);
            }
        }


        if (Pstream::master())
        {
            if (!obr_.foundObject<volScalarField>(pName_))
                return false;

            const volScalarField& pressure = obr_.lookupObject<volScalarField>(pName_);
            const fvMesh& mesh = pressure.mesh();
            dictionary forcesDict = dict_.subDict("forces");
            forcesDict.remove("CofR");
            forcesDict.add("CofR",vector::zero);

            functionObjects::forces f
            (
                "forces",
                mesh,
                forcesDict
            );

            f.calcForcesMoment();
            vector F = f.forceEff();  // fluid force
            vector M = f.momentEff();  // fluid moment about CofG

            outputVectorList[Pstream::myProcNo()].resize(faceCenters_.size(), vector::zero);
            vectorField toRot = faceCenters_;
            // structural velocity at fluid faces
            toRot = exampleModel_.currentFaceVelocity()();

            // quaternion R(0.0, 0.0, 0.0);
            quaternion R(quaternion::XYZ, vector::zero);
            // rotate velocities into global coordinates
            forAll(toRot, ii)
            {
                toRot[ii] = R.transform(toRot[ii]);
            }
            // send velocities in global coordinates, including rigid body forced velocity
            outputVectorList[Pstream::myProcNo()] = toRot;
        }
        else
        {
            outputVectorList[Pstream::myProcNo()].resize(faceCenters_.size(), vector::zero);

            // velocity at fluid faces
            vectorField toRot = exampleModel_.currentFaceVelocity()();

            // quaternion R(0.0, 0.0, 0.0);
            quaternion R(quaternion::XYZ, vector::zero);
            // rotate velocities into global coordinates
            forAll(toRot, ii)
            {
                toRot[ii] = R.transform(toRot[ii]);
            }
            // send velocity
            outputVectorList[Pstream::myProcNo()] = toRot;
        }
        Pstream::gatherList(outputVectorList);
    }
    else
    {
        outputVectorList[Pstream::myProcNo()].resize(faceCenters_.size(), vector::zero);
        Pstream::gatherList(outputVectorList);
    }

    // scatter
    vectorField outField(0,vector::zero);
    if (Pstream::master())
    {
        vectorField outputVectors(faceCenters_.size(), vector::zero);
        for (int ii = 0; ii < nExampleProc_; ii++)
            outputVectors += outputVectorList[ii];

        // scatter output field to different processors
        for (int ii = 0; ii < Pstream::nProcs(); ii++)
        {
            vectorField procOut(faceStarts_[ii+1]-faceStarts_[ii], vector::zero);
            for (int jj = 0; jj < procOut.size(); jj++)
                procOut[jj] = outputVectors[faceStarts_[ii]+jj];

            if (ii == Pstream::myProcNo())
                outField = procOut;
            else
            {
                UOPstream toNeighbor(ii, buf);
                toNeighbor << procOut;
            }
        }

        buf.finishedSends();
    }
    else
    {
        buf.finishedSends();
        UIPstream fromMaster(Pstream::masterNo(), buf);
        fromMaster >> outField;
    }

    forAll(patchIDs_, ii)
    {
        fvPatchField<vector>& patchField =
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>("U")
        ).boundaryField()[patchIDs_[ii]];

        List<vector> patchOutField(patchFaceStarts_[ii+1]-patchFaceStarts_[ii], vector::zero);
        for (int jj = patchFaceStarts_[ii]; jj < patchFaceStarts_[ii+1]; jj++)
            patchOutField[jj-patchFaceStarts_[ii]] = outField[jj];

        if (patchField.size() != patchOutField.size())
        {
            FatalErrorIn
            (
                "advanceFunctionObject::execute(const bool)"
            )   << "patchField and outField are not same size (" << patchField.size() << " versus " << patchOutField.size() << ") on processor " << Pstream::myProcNo()
                << exit(FatalError);
        }

        patchField.evaluate();
        patchField = patchOutField;
        patchField.updateCoeffs();
    }

    return true;
}


bool Foam::advanceFunctionObject::end()
{
    return true;
}


bool Foam::advanceFunctionObject::read(const dictionary& dict)
{
    Info << "MyProcNo = " << Pstream::myProcNo() << " of " << Pstream::nProcs() << endl;
    nExampleProc_ = Pstream::nProcs();
    return true;
}


// ************************************************************************* //
