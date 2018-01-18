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

#include "fvMeshFunctionObject.H"
#include "exampleClassModel.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class advanceFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class advanceFunctionObject
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Copy of the input dictionary
        dictionary dict_;

        //- Example control dictionary
        IOdictionary exampleDict_;

        //- Name of patch boundaries
        labelList patchIDs_;

        //- Switch for whether or not the time-step is adjustable
        bool adjustTimeStep_;

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

        //- Initialized or not. 1 = initialized, -1 attempted, but no "p" field
        int init_;


    // Protected Member Functions

        //- Return a description (type + name) for the output
        virtual string description() const;

        //- Get the dictionary to initialize the codeContext
        virtual const dictionary& exampleDict() const;


        //- Initialize
        bool startup();


        //- Disallow default bitwise copy construct
        advanceFunctionObject(const advanceFunctionObject&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const advanceFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("advance");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        advanceFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~advanceFunctionObject() = default;


    // Member Functions

        //- Read
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Write. Not used
        virtual bool write() { return true; }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
