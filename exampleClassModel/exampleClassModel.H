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
    Foam::exampleClassModel

Description

SourceFiles
    exampleModel.C

\*---------------------------------------------------------------------------*/

#ifndef exampleClassModel_H
#define exampleClassModel_H

#include "scalar.H"
#include "scalarMatrices.H"
#include "fvPatchField.H"
#include "pointField.H"
#include "vectorField.H"
#include "faceList.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class exampleClassModel;
Istream& operator>>(Istream&, exampleClassModel&);
Ostream& operator<<(Ostream&, const exampleClassModel&);


/*---------------------------------------------------------------------------*\
                      Class exampleClassModel Declaration
\*---------------------------------------------------------------------------*/

class exampleClassModel
{
    // Private data

        label nModes_;
        bool uncoupled_;
        labelListList fluidFaces2ExpMap_;
        List<List<vector> > fluidFaces2ExpWeights_;

        //- File path for time object
        fileName timePath_;

        //- Name of model
        word name_;

    // Private Member Functions

        int readModelData();

        tmp<Field<scalarField>> mapFluid2Example
        (
            const vectorField&,
            const vectorField&,
            vectorField&
        );

public:

    // Constructors

        exampleClassModel
        (
            word name,
            bool uncoupled,
            scalar& dt,
            fileName timePath
        );

    //- Destructor
    ~exampleClassModel() = default;


    // Member Functions

        //- Supply fluid points and faces to map
        tmp< Field<scalarField> > map2fluidDomain
        (
            const vectorField&,
            const vectorField&,
            vectorField&
        );

        //- Update the state vector for the iteration
        void updateState(scalarField&, vector, vector);

        //- Update the state vector for the iteration using input delta t
        void updateState(scalarField&, vector, vector, scalar);

        //- Update the state vector for the iteration and write output
        void updateState(scalarField&, vector, vector, word);

        //- Update the state vector for the iteration using input delta t and write output
        void updateState(scalarField&, vector, vector, scalar, word);

        //- Advace the state vectors for new time-step
        void advanceState(scalarField&, vector, vector);

        //- Advace the state vectors for new time-step using input delta t
        void advanceState(scalarField&, vector, vector, scalar);

        // Access

            //- Return the current velocities at the fluid faces
            tmp<vectorField> currentFaceVelocity();


        // IOstream Operators

        friend Istream& operator>>(Istream&, exampleClassModel&);
        friend Ostream& operator<<(Ostream&, const exampleClassModel&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
