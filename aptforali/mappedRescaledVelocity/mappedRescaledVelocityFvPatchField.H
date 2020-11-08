/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::mappedRescaledVelocityFvPatchField

Description
    Recycles the velocity and flux at a patch to this patch

    mode = NEARESTPATCHFACE : sample nearest face on selected patch
    mode = NEARESTFACE : sample nearest face on any patch. Note: does not
                         warn if nearest actually is on internal face!

SourceFiles
    mappedRescaledVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedRescaledVelocityFvPatchField_H
#define mappedRescaledVelocityFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "mappedFvPatch.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class mappedRescaledVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedRescaledVelocityFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- If true adjust the mapped field to maintain average value average_
        const bool setAverage_;

        //- Average value the mapped field is adjusted to maintain if
        //  setAverage_ is set true
        const vector average_;

        //- Interpolation scheme to use for nearestcell mode
        word interpolationScheme_;

        mutable autoPtr<interpolation<vector> > interpolator_;


        //- Boundary layer momentum thickness
        scalar theta_;

        //- Number of cells in wall-normal direction
        scalar nYCells_;

        //- Number of cells in spanwise direction
        scalar nZCells_;

        //- Max(nYCells, nZCells) //sadbg: redundant info ?? 
        scalar nCellsMax_;

        //- Y coordinate at the top wall, inflow boundary
        scalar YtopWall_;

        //- Y coordinate at the bottom wall, inflow boundary
        scalar YbottomWall_;

        //- Tolerance for indexing the faces
        scalar toleranceForIndexing_;

        //- Timescale to exclude transients
        scalar timeScaleToExcludeTransients_;

        //- Kinematic viscosity
        scalar kinematicViscosity_;

        //- delta99 BL thickness at inlet
        scalar delta99_;
    // Protected Member Functions

        //- Field to sample. Either on my or nbr mesh
        const volVectorField& sampleField() const;

        //- Access the interpolation method
        const interpolation<vector>& interpolator() const;

public:

    //- Runtime type information
    TypeName("mappedRescaledVelocityNew");


    // Constructors

        //- Construct from patch and internal field
        mappedRescaledVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedRescaledVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mappedRescaledVelocityFvPatchField
        //  onto a new patch
        mappedRescaledVelocityFvPatchField
        (
            const mappedRescaledVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedRescaledVelocityFvPatchField
        (
            const mappedRescaledVelocityFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new mappedRescaledVelocityFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedRescaledVelocityFvPatchField
        (
            const mappedRescaledVelocityFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new mappedRescaledVelocityFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
