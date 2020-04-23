/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "kEpsilonCC.H"
#include "fvOptions.H"
#include "bound.H"
//#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilonCC<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonCC<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonCC<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonCC<BasicTurbulenceModel>::kEpsilonCC
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    // Original constants
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

// curvature correction
    cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr1",
            this->coeffDict_,
            1.0
        )
    ),
    cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr2",
            this->coeffDict_,
            2.0
        )
    ),
    cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr3",
            this->coeffDict_,
            1.0
        )
    ),
// end curvature correction

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )

{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonCC<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        // Curvature correction terms
        cr1_.readIfPresent(this->coeffDict());
        cr2_.readIfPresent(this->coeffDict());
        cr3_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kEpsilonCC<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    // Algebraic terms 
    tmp<volTensorField> tgradU = fvc::grad(U);
    tmp<volTensorField> tSkew = skew(tgradU()); 
    tmp<volSymmTensorField> tSymm = symm(tgradU());  
    volScalarField S2(2*magSqr(tSymm()));        
//
    //volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    //volScalarField::Internal G(this->GName(), nut()*GbyNu);
    //tgradU.clear();

    // Curvature Correction terms
    // Determination of rStar
    // Strain Rate Tensor S
    volScalarField symInnerProduct(2.0*tSymm() && tSymm()); 
    // Vorticity tensor Omega
    volScalarField asymInnerProduct
    (
        max(2.0*tSkew() && tSkew(),
        dimensionedScalar("0", dimensionSet(0, 0, -2, 0, 0), 0.0))
    );

    volScalarField w
    (
        atan(dimensionedScalar("4",dimensionSet(0,0,2,0,0),1.0e-02)*asymInnerProduct)*2.0/(constant::mathematical::pi)*(asymInnerProduct-symInnerProduct)+symInnerProduct
    );
    // rStar
    // Avoid Division by zero with SMALL keymord (?)
    volScalarField rStar(sqrt(symInnerProduct/max(w, dimensionedScalar("minw",
                                                                       w.dimensions(),
                                                                       SMALL))));
    // volScalarField rStar(sqrt(symInnerProduct/w));
    // Determination of rTilda
    // using omega = epsilon/k to replace omega
    volScalarField D(sqrt(max(symInnerProduct, epsilon_*epsilon_/k_/k_/0.09)));
    tmp<volSymmTensorField> divS =
    (
        fvc::ddt(tSymm())
       +fvc::div
        (
            alphaRhoPhi, tSymm()
        )
    );
    volScalarField rT((tSkew().T() & tSymm) && divS);

    divS.clear();
    tSkew.clear();
    tSymm.clear();
    
    volScalarField w2
    (
        atan(dimensionedScalar("1",dimensionSet(0,0,2,0,0),1.0e-2)*asymInnerProduct)*2.0/(constant::mathematical::pi)*(sqrt(asymInnerProduct)-D)+D //T
    );
    volScalarField rTilda(2.0*rT/w2/D/D/D);
   
    // ////
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    volScalarField Fr1    //curvature correction
    (
	max
    	(
		min
		(
        		(scalar(1.0)+cr1_)*2.0*rStar/(scalar(1)+rStar)*(scalar(1.0)
			-cr3_*atan(cr2_*rTilda))-cr1_,
        		scalar(1.25)
		),
			scalar(0.0)
	)
    );

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()*Fr1()/k_() // added Fr1 to correct curvature
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G*Fr1() // curvature correction
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
