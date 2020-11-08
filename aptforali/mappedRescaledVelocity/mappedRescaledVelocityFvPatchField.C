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

\*---------------------------------------------------------------------------*/

#include "createFields.H"
#include "mappedRescaledVelocityFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedRescaledVelocityFvPatchField::
mappedRescaledVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    setAverage_(false),
    average_(pTraits<vector>::zero),
    interpolationScheme_(interpolationCell<vector>::typeName),
    theta_(0.0072),
    nYCells_(95),
    nZCells_(95),
    nCellsMax_(95),
    toleranceForIndexing_(1.e-07),
    kinematicViscosity_(3.697e-05),
    delta99_(0.073)
{}


Foam::mappedRescaledVelocityFvPatchField::
mappedRescaledVelocityFvPatchField
(
    const mappedRescaledVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_),
    interpolationScheme_(ptf.interpolationScheme_),
    theta_(ptf.theta_),
    nYCells_(ptf.nYCells_),
    nZCells_(ptf.nZCells_),
    nCellsMax_(ptf.nCellsMax_),
    toleranceForIndexing_(ptf.toleranceForIndexing_),
    kinematicViscosity_(ptf.kinematicViscosity_),
    delta99_(ptf.delta99_)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "mappedRescaledVelocityFvPatchField::"
            "mappedRescaledVelocityFvPatchField"
            "("
                "const mappedRescaledVelocityFvPatchField&, "
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "const fvPatchFieldMapper&"
            ")"
        )   << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


Foam::mappedRescaledVelocityFvPatchField::
mappedRescaledVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    setAverage_(readBool(dict.lookup("setAverage"))),
    average_(pTraits<vector>(dict.lookup("average"))),
    interpolationScheme_(interpolationCell<vector>::typeName),
    theta_(readScalar(dict.lookup("momentumThickness"))),
    nYCells_(readScalar(dict.lookup("nCellsWallNormal"))),
    nZCells_(readScalar(dict.lookup("nCellsSpanwise"))),
    nCellsMax_(readScalar(dict.lookup("nCellsMax"))),
    toleranceForIndexing_(readScalar(dict.lookup("toleranceForIndexing"))),
    kinematicViscosity_(readScalar(dict.lookup("kinematicViscosity"))),
    delta99_(readScalar(dict.lookup("delta99")))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "mappedRescaledVelocityFvPatchField::"
            "mappedRescaledVelocityFvPatchField"
            "("
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "const dictionary&"
            ")"
        )   << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        this->patch().patch()
    );
    if (mpp.mode() == mappedPatchBase::NEARESTCELL)
    {
        dict.lookup("interpolationScheme") >> interpolationScheme_;
    }

    // if (mpp.mode() == mappedPolyPatch::NEARESTCELL)
    // {
    //     FatalErrorIn
    //     (
    //         "mappedRescaledVelocityFvPatchField::"
    //         "mappedRescaledVelocityFvPatchField"
    //         "("
    //             "const fvPatch&, "
    //             "const DimensionedField<vector, volMesh>&, "
    //             "const dictionary&"
    //         ")"
    //     )   << "Patch " << p.name()
    //         << " of type '" << p.type()
    //         << "' can not be used in 'nearestCell' mode"
    //         << " of field " << dimensionedInternalField().name()
    //         << " in file " << dimensionedInternalField().objectPath()
    //         << exit(FatalError);
    // }
}


Foam::mappedRescaledVelocityFvPatchField::
mappedRescaledVelocityFvPatchField
(
    const mappedRescaledVelocityFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_),
    interpolationScheme_(ptf.interpolationScheme_),
    theta_(ptf.theta_),
    nYCells_(ptf.nYCells_),
    nZCells_(ptf.nZCells_),
    nCellsMax_(ptf.nCellsMax_),
    toleranceForIndexing_(ptf.toleranceForIndexing_),
    kinematicViscosity_(ptf.kinematicViscosity_),
    delta99_(ptf.delta99_)
{}


Foam::mappedRescaledVelocityFvPatchField::
mappedRescaledVelocityFvPatchField
(
    const mappedRescaledVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_),
    interpolationScheme_(ptf.interpolationScheme_),
    theta_(ptf.theta_),
    nYCells_(ptf.nYCells_),
    nZCells_(ptf.nZCells_),
    nCellsMax_(ptf.nCellsMax_),
    toleranceForIndexing_(ptf.toleranceForIndexing_),
    kinematicViscosity_(ptf.kinematicViscosity_),
    delta99_(ptf.delta99_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::mappedRescaledVelocityFvPatchField::sampleField() const
{

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
     mappedRescaledVelocityFvPatchField::patch().patch()
    );
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const word& fieldName = dimensionedInternalField().name();
    
    return nbrMesh.lookupObject<volVectorField>(fieldName);
}

const Foam::interpolation<Foam::vector>& Foam::mappedRescaledVelocityFvPatchField::interpolator() const
{
    if (!interpolator_.valid())
    {
        interpolator_ = interpolation<vector>::New
        (
            interpolationScheme_,
            sampleField()
        );
    }
    return interpolator_();
}

void Foam::mappedRescaledVelocityFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
     mappedRescaledVelocityFvPatchField::patch().patch()
    );
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const word& fieldName = dimensionedInternalField().name();
    const volVectorField& UField =
        nbrMesh.lookupObject<volVectorField>(fieldName);

    vectorField newValues;

    switch (mpp.mode())
    {

        case mappedPatchBase::NEARESTCELL:
        {
            const mapDistribute& distMap = mpp.map();

            if (interpolationScheme_ != interpolationCell<vector>::typeName)
            {
                // Send back sample points to the processor that holds the cell
                vectorField samples(mpp.samplePoints());
                distMap.reverseDistribute
                (
		 (mpp.sameRegion() ? this->patch().boundaryMesh().mesh().nCells() : nbrMesh.nCells()),
                    point::max,
                    samples
                );

                const interpolation<vector>& interp = interpolator();

                newValues.setSize(samples.size(), pTraits<vector>::max);
                forAll(samples, cellI)
                {
                    if (samples[cellI] != point::max)
                    {
                        newValues[cellI] = interp.interpolate
                        (
                            samples[cellI],
                            cellI
                        );
                    }
                }
            }
            else
            {
                newValues = sampleField();
                //cout << "I am here ... " << "\n";
            }

            //cout << "I am here ... " << "\t" << newValues.size()<< "\n";
            distMap.distribute(newValues);
            //cout << "I am here ... " << "\t" << newValues.size()<< "\n";

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "mappedRescaledVelocityFvPatchField::"
                "updateCoeffs()"
            )   << "patch can only be used in NEARESTCELL, "
                << "mode" << nl
                << abort(FatalError);
        }
    }

     // * * * * * * * * * * * * * * * sadbg: Rescaling the velocity data  * * * * * * * * * * * * * //
     //sadbg: if (rescale_)
     //sadbg: Take inlet boundary layer thickness targetted as an user input and may be nmax as well
     //sadbg: Note that the spanwise averaging assumes sequential storage of cells and needs to be tested for curvilinear geometries before using this
     //sadbg: Take freestream velocity as input
    int rank                        =  Pstream::myProcNo();

//    const volVectorField& UField    =  nbrMesh.lookupObject<volVectorField>(fldName);
    vectorField Urecycled           =  UField.internalField();
    scalar deltaT_                  =  UField.time().deltaT().value();
    if(rank == 0)
      {
        t_  = UField.time().value();

        if(counter==1)
          {
            t0_ = UField.time().value();
          }

        T_  = (t_ - t0_)+1e-08;

        nSteps_ = (T_/deltaT_) + 1;

	int     i, j, m, n;
        //const int nmax = 103, nCellsWallNormal = 95, nCellsSpanwise = 103;
        const int nmax = nCellsMax_, nCellsWallNormal = nYCells_, nCellsSpanwise = nZCells_;
	
	//const int nmax = 95, nCellsWallNormal = 95, nCellsSpanwise = 45;
    
    scalar Ybottom = YbottomWall_; 

	scalar urecy[nmax][nmax], vrecy[nmax][nmax], wrecy[nmax][nmax];
	scalar uInstantInlet[nmax][nmax], vInstantInlet[nmax][nmax], wInstantInlet[nmax][nmax];
	
	m = nCellsSpanwise-1; //minus 1 because cell index starts from zero
	n = nCellsWallNormal-1; //minus 1 because cell index starts from zero

	scalar  Uaveraged[nmax], totVol[nmax];
	scalar  Uinfty;
	scalar  thetaTarget, thetaInlet;
	scalar ytemp, ztemp;
	scalar Yrecy[nmax];

	//sadbg: Get thetaInlet as input parameters
	thetaTarget = theta_; //0.000635; //0.635; 
	
	//Intialize
	thetaInlet = thetaTarget;

	//sadbg: Initialize Uaveraged
	for(i=0; i<nmax; i++)
	  {
	    Uaveraged[i] = 0.0;
	    totVol[i] = 0.0;
	  }

	label inletPatchID = nbrMesh.boundaryMesh().findPatchID("INLET");
	const vectorField& cc1 = nbrMesh.Cf().boundaryField()[inletPatchID];
	const scalarField& vol1 = nbrMesh.magSf().boundaryField()[inletPatchID];
	scalarField xx = cc1.component(vector::X);
	scalarField yy = cc1.component(vector::Y);
	scalarField zz = cc1.component(vector::Z);
	scalarField rr = yy; //zz;
	
	//sadbg: Build face list in a sorted order for radius (which is z-coordinate for this case)
	i = 0;
	j = 0;
	int count;
	bool indexed = false;
        scalar tol_ = toleranceForIndexing_; //1.e-07; //1.e-07; //1.e-03; //sadbg: TAKE THIS AS INPUT

        //cout << "I am here ... " <<"\t" << newValues.size() << "\n";

	forAll(newValues, faceI)
	  {
	    //sadbg: Check for already indexed faces
	    //sadbg: Take the tolerance as user input to generalize this code
	    for(count=0;count<i;count++)
	      {
		if(mag(storedRad[count]-rr[faceI]) <= tol_)
		  {
		    // cout << "I am here ... " <<"\t" << i << "\t" << rr[faceI] << "\n";
		    indexed = true;
		  }
	      }

	    if(indexed==false)
	      {
		j = 0;
		indexI[faceI] = i;
		indexJ[faceI] = j;
		if(faceI > 50000)
			cout<<"error is here"<<endl;
		forAll(newValues, faceJ)
		  {
		    if((faceI != faceJ)&&(mag(rr[faceJ]-rr[faceI]) <= tol_))
		      {
                        j++;
			indexI[faceJ] = i;
			indexJ[faceJ] = j;
			storedRad[i] = rr[faceJ];
				      
		}
		  }
		i++;
	      }

	    indexed = false;
	  }
	
        //sadbg: It is not clear if the mapped values are stored on the inflow patch
	// in the same way as sample plane. But, that is an assumption
	// Did not give any problem for annulus.
	forAll(newValues, faceI)
	  {
	    ytemp = yy[0];
	    ztemp = zz[0];
	    Yrecy[0] = ytemp; //ztemp;
	    
	    i = indexI[faceI];
	    j = indexJ[faceI];
	    Yrecy[i] = rr[faceI];
	    Urecycled[faceI] = newValues[faceI]; //sadbg: KEY STEP
	    
            //Info << "i = "<< i<< "j = "<< j<<"face Num = "<< faceI<<endl;

	    urecy[i][j]  = Urecycled[faceI].x(); 
	    vrecy[i][j]  = Urecycled[faceI].y();
	    wrecy[i][j]  = Urecycled[faceI].z();

	    Uaveraged[i] = Uaveraged[i] + urecy[i][j]*vol1[faceI];
	    totVol[i]    = totVol[i] + vol1[faceI];
	  
            //cout << i << "\t" << totVol[i] << "\n";
          }
	
	//std::ofstream profile("profile.dat", ios_base::app);
	//std::ofstream UaveragedProfile("UaveragedProfile.dat");
        //std::ofstream thetaWrite("thetaWrite.dat");

 	// sadbg: Average in the spanwise direction
 	for(i=0; i<=n; i++)
 	  {
 	    Uaveraged[i] = Uaveraged[i]/(totVol[i]+1.e-12);
	//Info << Uaveraged[i] << endl;
 	  }
	Uinfty = Uaveraged[n];
	//Info << "Uinfty = " << Uinfty << endl;
//
    //    for(i=0; i<=n; i++)
    //      {
    //        UaveragedProfile << i << "\t" << Yrecy[i] << "\t"<< Uaveraged[i] << "\n";
    //      }
	//
	//sadbg: Find momentum thickness
 	thetaRecyBottom = 0.0;
 	for(i=0; i<n; i++)
 	  {
 	    thetaRecyBottom = thetaRecyBottom + (Uinfty - Uaveraged[i])*(Uaveraged[i]/Uinfty)*(Yrecy[i+1]-Yrecy[i])/Uinfty;
 	  }

        //thetaWrite << thetaRecyBottom ;
     
	//sadbg: Composite formula to rescale the inflow velocity field stored in newValues
	scalar utemp[nmax][nmax], vtemp[nmax][nmax], wtemp[nmax][nmax];
	scalar zshiftU[nmax][nmax], zshiftV[nmax][nmax], zshiftW[nmax][nmax];
	
	for(j=0;j<=m;j++)
	  {
	    for(i=0;i<=n;i++)
	      {
		utemp[i][j] = urecy[i][j];
		vtemp[i][j] = vrecy[i][j];
		wtemp[i][j] = wrecy[i][j];
	      }
	  }
	
	//sadbg: Disorganize spanwise durable structures by mirroring
        for(j=0; j<=m; j++)
          {
	    for(i=0; i<=n; i++)
	      {
       		zshiftU[i][m-j] = utemp[i][j];
        	zshiftV[i][m-j] = vtemp[i][j];
       	 	zshiftW[i][m-j] = -wtemp[i][j];
	      }
	  }


	//sadbg: Compute Yplus at inlet and recycling stations
	scalar Yinlt_[nmax];
	    
	for(i=0;i<=n;i++)
	  {
	    Yinlt_[i]  = (Yrecy[i])*thetaRecyBottom/thetaInlet;
	  }

	for(i=0;i<=n;i++)
	  {
	    for(j=0;j<=m;j++)
	      {
		
		bool updated_ = false;

		if((Yinlt_[i]-Yrecy[n]) >= 1.e-08)
		  {
		    uInstantInlet[i][j] = zshiftU[n][j];
		    vInstantInlet[i][j] = zshiftV[n][j];
		    wInstantInlet[i][j] = zshiftW[n][j];
		    
		    updated_ = true;
		  }
		else
		  {
		    for(int k=0; k<(n-1); k++)
		      {
			if((updated_==false)&&((Yinlt_[i]-Yrecy[k]) > 1.e-08) && ((Yinlt_[i]-Yrecy[k+1]) < 1.e-08))
			  {
			    scalar L0, L1, L2;
			    
			    L0 = (Yinlt_[i]-Yrecy[k+1])*(Yinlt_[i]-Yrecy[k+2])/((Yrecy[k]-Yrecy[k+1])*(Yrecy[k]-Yrecy[k+2]));
			    L1 = (Yinlt_[i]-Yrecy[k])*(Yinlt_[i]-Yrecy[k+2])/((Yrecy[k+1]-Yrecy[k])*(Yrecy[k+1]-Yrecy[k+2]));
			    L2 = (Yinlt_[i]-Yrecy[k])*(Yinlt_[i]-Yrecy[k+1])/((Yrecy[k+2]-Yrecy[k])*(Yrecy[k+2]-Yrecy[k+1]));
			    
			    uInstantInlet[i][j] = zshiftU[k][j]*L0 + zshiftU[k+1][j]*L1 + zshiftU[k+2][j]*L2;
			    vInstantInlet[i][j] = zshiftV[k][j]*L0 + zshiftV[k+1][j]*L1 + zshiftV[k+2][j]*L2;
			    wInstantInlet[i][j] = zshiftW[k][j]*L0 + zshiftW[k+1][j]*L1 + zshiftW[k+2][j]*L2;

			    updated_ = true;
			    break;
			  }
		      }
		  }
		
		if(updated_ == false)
		  {
		    uInstantInlet[i][j] = zshiftU[i][j];
		    vInstantInlet[i][j] = zshiftV[i][j];
		    wInstantInlet[i][j] = zshiftW[i][j];
		  }
		
	      }
	  }

	//sadbg: Assign scaled values
	forAll(newValues, faceI)
	  {
	    i = indexI[faceI];
	    j = indexJ[faceI];


	    newValues[faceI].x() = uInstantInlet[i][j];
	    newValues[faceI].y() = vInstantInlet[i][j];
	    newValues[faceI].z() = wInstantInlet[i][j];
	
       if(Yrecy[i] > (Ybottom + delta99_))
           {
                newValues[faceI].x() = Uinfty;
                newValues[faceI].y() = 0.0;
                newValues[faceI].z() = 0.0;
            }
	
	  }
	Info << "Uinfty = " << Uinfty << endl;
	//profile << t_ << "\t" << Uinfty * thetaRecyBottom /(1.668e-5) << "\t" << thetaInlet << "\n";

        tcheck_ = t_;
        counter++;
      }

    // * * * * * * * * * * * * * * * sadbg: Rescaling the velocity data  * * * * * * * * * * * * * //

    if (setAverage_)
    {
        Foam::vector averagePsi =
            gSum(this->patch().magSf()*newValues)
           /gSum(this->patch().magSf());

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            newValues *= mag(average_)/mag(averagePsi);
        }
        else
        {
            newValues += (average_ - averagePsi);
        }
    }

    operator==(newValues);

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::mappedRescaledVelocityFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
    os.writeKeyword("average") << average_ << token::END_STATEMENT << nl;
    os.writeKeyword("interpolationScheme") << interpolationScheme_
        << token::END_STATEMENT << nl;
    os.writeKeyword("momentumThickness") << theta_ << token::END_STATEMENT << nl;
    os.writeKeyword("nCellsWallNormal") << nYCells_ << token::END_STATEMENT << nl;
    os.writeKeyword("nCellsSpanwise") << nZCells_ << token::END_STATEMENT << nl;
    os.writeKeyword("nCellsMax") << nCellsMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("toleranceForIndexing") << toleranceForIndexing_ << token::END_STATEMENT << nl;
    os.writeKeyword("kinematicViscosity") << kinematicViscosity_ << token::END_STATEMENT << nl;
    os.writeKeyword("delta99") << delta99_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        mappedRescaledVelocityFvPatchField
    );
}


// ************************************************************************* //
