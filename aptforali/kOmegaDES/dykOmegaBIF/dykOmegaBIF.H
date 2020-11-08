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
    Foam::incompressible::RASModels::dykOmegaBIF

Implementation for 
Yin, Zifei, and Paul A. Durbin. "An adaptive DES smodel that allows 
wall-resolved eddy simulation." 
International Journal of Heat and Fluid Flow (2016).

Any questions, please contact Zifei Yin at zifeiyin@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef dykOmegaBIF_H
#define dykOmegaBIF_H

#include "RASModel.H"
#include "wallDist.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class dykOmegaBIF Declaration
\*---------------------------------------------------------------------------*/

class dykOmegaBIF
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar beta_;
            dimensionedScalar alpha_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;
            dimensionedScalar Cu_;
	    dimensionedScalar kappa_;
            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;


            dimensionedScalar betaStar_;

            dimensionedScalar alpha1_;
            dimensionedScalar alpha2_;
            dimensionedScalar alpha3_;
            dimensionedScalar alpha4_;
            dimensionedScalar alpha5_;
            dimensionedScalar Cr_;
            dimensionedTensor rotaxis_;
            dimensionedScalar rotspeed_;//rad/s
            dimensionedScalar unit1_;
            dimensionedScalar unit2_;
            dimensionedScalar unit3_;
            dimensionedScalar unit4_;
            dimensionedScalar m_;


        // Fields
            volSymmTensorField Sij_;
            volScalarField I1_;
            volScalarField I2_;
            volScalarField I3_;
            volScalarField IIS_;
            volScalarField IIIS_;
            volScalarField k_;
            volScalarField omega_;
            volScalarField nut_;
	    volScalarField lDDES_;
	    volScalarField lRANS_;
	    volScalarField lLES_;
	    volScalarField lzero_;
	    volScalarField fd_;
	    volScalarField rd_;

	    volScalarField delta_;
	    volScalarField cube_root_vol_;
	    volScalarField hmax_;

            volScalarField CDES_;

        volScalarField Lk_;
        volScalarField Fr_;
        volScalarField Cbound_;
        volScalarField epsilon_tmp;
	    wallDist y_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;

	volScalarField cD(const volSymmTensorField& D, const volScalarField& w) const;
    // Protected Member Functions
        tmp<volScalarField> f() const;
public:

    //- Runtime type information
    TypeName("dykOmegaBIF");

    // Constructors

        //- Construct from components
        dykOmegaBIF
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~dykOmegaBIF()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff",(alphaK_*k_/omega_ + nu()))
//                new volScalarField("DkEff", alphaK_*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff",  (alphaOmega_*k_/omega_ + nu()))
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return delta
        virtual tmp<volScalarField> delta() const
        {
            return delta_;
        }

	// Calculate delta
	void calcdelta();

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();

	void update_nut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
