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

Class
    Foam::advanceFunctionObject

Description
    functionObject using dynamic code compilation.

SourceFiles
    advanceFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef advanceFunctionObject_H
#define advanceFunctionObject_H

#include "functionObject.H"
#include "exampleClassModel.H"
#include "Switch.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class advanceFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class advanceFunctionObject
:
    public functionObject
{
protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Input dictionary
        dictionary dict_;

        //- Example control dictionary
        IOdictionary exampleDict_;

        //- Object registry containing all fields
        const objectRegistry& obr_;

        //- Name of patch boundaries
        labelList patchIDs_;

        //- Switch for whether or not the time-step is adjustable
        Switch adjustTimeStep_;

        //- Name of the pressure field
        word pName_;

        //- Underlying example class model object
        exampleClassModel exampleModel_;

        //- Start of faces for each processor in agglomerated list. Only used if master process
        labelList faceStarts_;

        //- Start of points for each processor in agglomerated list. Only used if master process
        labelList pointStarts_;

        //- Start of faces for each patch in agglomerated processor list. Unique for each processor
        labelList patchFaceStarts_;

        //- Start of points for each patch in agglomerated processor list. Unique for each processor
        labelList patchPointStarts_;

        //- Face centers for entire combined patch - used for approximate BC
        vectorField faceCenters_;

        //- Points for entire combined patch - used for exact BC
        vectorField fluidPts_;

        //- Output interval - # of timesteps
        label outputInterval_;

        //- Current time index
        label curTimeIndex_;

        //- Number of processors to split between
        label nExampleProc_;

    // Protected Member Functions

        //- Return a description (type + name) for the output
        virtual string description() const;

        //- Get the dictionary to initialize the codeContext
        virtual const dictionary& exampleDict() const;

private:

        //- Disallow default bitwise copy construct
        advanceFunctionObject(const advanceFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const advanceFunctionObject&);


public:

    //- Runtime type information
    TypeName("advance");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        advanceFunctionObject
        (
            const word& name,
            const Time& time,
            const dictionary& dict,
            bool readNow=true   // allow child-classes to avoid compilation
        );


    //- Destructor
    virtual ~advanceFunctionObject();

    // Member Functions

        //- Called at the start of the time-loop
        virtual bool start();

        //- Called at each ++ or += of the time-loop. forceWrite overrides the outputControl behaviour.
        virtual bool execute(const bool forceWrite);

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);

        virtual void updateMesh(const Foam::mapPolyMesh&) {}

        virtual void movePoints(const Foam::polyMesh&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
