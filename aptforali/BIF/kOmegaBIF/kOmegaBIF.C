/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "kOmegaBIF.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaBIF, 0);
addToRunTimeSelectionTable(RASModel, kOmegaBIF, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kOmegaBIF::f() const
{
     volSymmTensorField Sij_ (symm(fvc::grad(U_)));

     /*volScalarField IIS = Sij_&&Sij_;//s^-2
     volScalarField IIIS = Sij_&Sij_&&Sij_;//s^-3*/
volScalarField I1 = tr(Sij_);//s^_1
     volScalarField I2 = Sij_&&Sij_;//s^-2
     volScalarField I3 = Sij_&Sij_&&Sij_;//s^-3
     volScalarField IIS = 0.5*(sqr(I1)-I2);//s^-2
     volScalarField IIIS = (pow(I1,3)-3*I1*I2+2*I3)/6;//s^-3

     volSymmTensorField Xij = (Sij_-Sij_)*unit4_;

Xij.replace(0,(sqr(IIS)+12*IIIS*Sij_.component(0)+6*IIS*(Sij_.component(0)*Sij_.component(0)+Sij_.component(1)*Sij_.component(1)+Sij_.component(2)*Sij_.component(2)))/max(2*pow(IIS,3)-12*sqr(IIIS),pow(min(max(k_,kMin_),kMin_)/unit2_,3)));

Xij.replace(1,(12*IIIS*Sij_.component(1)+6*IIS*(Sij_.component(0)*Sij_.component(1)+Sij_.component(1)*Sij_.component(3)+Sij_.component(2)*Sij_.component(4)))/max(2*pow(IIS,3)-12*sqr(IIIS),pow(min(max(k_,kMin_),kMin_)/unit2_,3)));
     
Xij.replace(2,(12*IIIS*Sij_.component(2)+6*IIS*(Sij_.component(0)*Sij_.component(2)+Sij_.component(1)*Sij_.component(4)+Sij_.component(2)*Sij_.component(5)))/max(2*pow(IIS,3)-12*sqr(IIIS),pow(min(max(k_,kMin_),kMin_)/unit2_,3)));

Xij.replace(3,(sqr(IIS)+12*IIIS*Sij_.component(3)+6*IIS*(Sij_.component(1)*Sij_.component(1)+Sij_.component(3)*Sij_.component(3)+Sij_.component(4)*Sij_.component(4)))/max(2*pow(IIS,3)-12*sqr(IIIS),pow(min(max(k_,kMin_),kMin_)/unit2_,3)));

Xij.replace(4,(12*IIIS*Sij_.component(4)+6*IIS*(Sij_.component(1)*Sij_.component(2)+Sij_.component(3)*Sij_.component(4)+Sij_.component(4)*Sij_.component(5)))/max(2*pow(IIS,3)-12*sqr(IIIS),pow(min(max(k_,kMin_),kMin_)/unit2_,3)));

Xij.replace(5,(sqr(IIS)+12*IIIS*Sij_.component(5)+6*IIS*(Sij_.component(2)*Sij_.component(2)+Sij_.component(4)*Sij_.component(4)+Sij_.component(5)*Sij_.component(5)))/max(2*pow(IIS,3)-12*sqr(IIIS),pow(min(max(k_,kMin_),kMin_)/unit2_,3)));


volSymmTensorField divS = fvc::ddt(Sij_) + fvc::div(surfaceScalarField("phiU",phi_), Sij_);
    volTensorField OmegaF = rotaxis_*rotspeed_*unit1_*k_/k_;
   volTensorField OmegaSS = OmegaF-((Sij_&divS)-(divS&Sij_))/max(2*sqr(mag(Sij_)),unit3_);

    volVectorField FeSe = (U_-U_)/y_;
    FeSe.replace(0, OmegaF.component(5)-OmegaF.component(7)+ OmegaSS.component(7)-OmegaSS.component(5)); 
    FeSe.replace(1, OmegaF.component(6)-OmegaF.component(2)+ OmegaSS.component(2)-OmegaSS.component(6)); 
    FeSe.replace(2, OmegaF.component(1)-OmegaF.component(3)+ OmegaSS.component(3)-OmegaSS.component(1));   

    volVectorField w = (U_-U_)/y_;
    w.replace(0,IIS*(Xij.component(0)*FeSe.component(0)+Xij.component(1)*FeSe.component(1)+Xij.component(2)*FeSe.component(2)));
    w.replace(1,IIS*(Xij.component(1)*FeSe.component(0)+Xij.component(3)*FeSe.component(1)+Xij.component(4)*FeSe.component(2)));
    w.replace(2,IIS*(Xij.component(2)*FeSe.component(0)+Xij.component(4)*FeSe.component(1)+Xij.component(5)*FeSe.component(2)));

    volTensorField neijkwi = (-fvc::grad(U_)+fvc::grad(U_));
    neijkwi.replace(1,-w.component(2));
    neijkwi.replace(2,w.component(1));
    neijkwi.replace(3,w.component(2));
    neijkwi.replace(5,-w.component(0));
    neijkwi.replace(6,-w.component(1));
    neijkwi.replace(7,w.component(0));

 
     volTensorField Wij = OmegaF+neijkwi;

    volTensorField Omega_mod = -skew(fvc::grad(U_))+OmegaF+(Cr_-1)*Wij;

    volScalarField T1 = scalar(1)/(betaStar_*omega_);
    volScalarField T2 = scalar(6)*sqrt(nu()/(betaStar_*omega_*k_));
    volScalarField T3 = pow(pow(T1,1.625)*T2,1/2.625);
    volScalarField Tm = max(T1,T3);  
    volScalarField yita1 = magSqr(Sij_)*sqr(Tm);
    volScalarField yita2 = magSqr(Omega_mod)*sqr(Tm);
    volScalarField yita3 = yita1-yita2; 
    //volScalarField f = 1/(alpha1_*(mag(yita3)-yita3)+sqrt(1-min(alpha2_*yita3,0.99)));//correction factor for the Cmu_  

	return min(betaStar_*max(pow(sqrt((1+alpha5_*yita1)/(1+alpha5_*yita2))+alpha1_*sqrt(yita2)*sqrt(mag(yita3)-yita3),-1),0.000001)*(1+alpha2_*mag(yita3)+alpha3_*yita3)/(1+alpha4_*mag(yita3)),2.5);

   // return min(1/(alpha1_*(mag(yita3)-yita3)+sqrt(1-min(alpha2_*yita3,scalar(0.99)))),2.5);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaBIF::kOmegaBIF
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            0.072
        )
    ),
    alpha_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha",
            coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.5
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    alpha1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha1",
            coeffDict_,
            0.062
        )
    ),
    alpha2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha2",
            coeffDict_,
            0.8
        )
    ),
    alpha3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha3",
            coeffDict_,
            0.4
        )
    ),
   alpha4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha4",
            coeffDict_,
            0.4
        )
    ),
   alpha5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha5",
            coeffDict_,
            0.025
        )
    ),

    Cr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr",
            coeffDict_,
            2
        )
    ),

        rotaxis_
    (
        dimensioned<tensor>::lookupOrAddToDict
        (
            "rotaxis",
            coeffDict_,
            tensor(0,-1,0,1,0,0,0,0,0)
        )
    ),
    rotspeed_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "rotspeed",
            coeffDict_,
            0
        )
    ),
    unit1_
    (
        dimensioned<scalar>
        (
            "unit1_",
             dimensionSet(0,0,-1,0,0,0,0),
             1
        )
    ),
    unit2_
    (
        dimensioned<scalar>
        (
            "unit2_",
             dimensionSet(0,2,0,0,0,0,0),
             1
        )
    ),
    unit3_
    (
        dimensioned<scalar>
        (
            "unit3_",
             dimensionSet(0,0,-2,0,0,0,0),
             1e-300
        )
    ),
    unit4_
    (
        dimensioned<scalar>
        (
            "unit4_",
             dimensionSet(0,0,3,0,0,0,0),
             1
        )
    ),
    m_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "m",
            coeffDict_,
            3
        )
    ),

    y_(mesh_),
	Sij_ (symm(fvc::grad(U_))),
	I1_ (tr(Sij_)),
	I2_ (Sij_&&Sij_),
	I3_ (Sij_&Sij_&&Sij_),
	IIS_ (0.5*(sqr(I1_)-I2_)),
	IIIS_ ((pow(I1_,3)-3*I1_*I2_+2*I3_)/6),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    nut_ = f()*k_/omega_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaBIF::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kOmegaBIF::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kOmegaBIF::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> kOmegaBIF::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool kOmegaBIF::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());
        alphaK_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());
        rotaxis_.readIfPresent(coeffDict());
        rotspeed_.readIfPresent(coeffDict());
        m_.readIfPresent(coeffDict());
        alpha1_.readIfPresent(coeffDict());
        alpha2_.readIfPresent(coeffDict());
        alpha3_.readIfPresent(coeffDict());
        alpha4_.readIfPresent(coeffDict());
        alpha5_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaBIF::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));
   	
	volScalarField fr
	(
        IOobject
	(
            "fr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        f()
    );
    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
        alpha_*G*omega_/k_
      - fvm::Sp(beta_*omega_, omega_)
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(Cmu_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = f()*k_/omega_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
