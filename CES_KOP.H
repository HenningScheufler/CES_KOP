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
    Foam::incompressible::HybridLesRasModels::CES_KOS

Description
    Linear unifified model with new damping Function

SourceFiles
    CES_KOS.C

\*---------------------------------------------------------------------------*/

#ifndef CES_KOS_H
#define CES_KOS_H

#include "HybridModel.H"
#include "LESfilter.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace HybridLesRasModels
{

/*---------------------------------------------------------------------------*\
                        Class CES_KOS Declaration
\*---------------------------------------------------------------------------*/

class CES_KOS
:
    public HybridLesRasModel
{
    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields();

        // Disallow default bitwise copy construct and assignment
        CES_KOS(const CES_KOS&);
        CES_KOS& operator=(const CES_KOS&);


protected:

    // Protected data

        // Model constants
	dimensionedScalar cOmega1_;
	dimensionedScalar cOmega2_;
	dimensionedScalar cOmega_;
	dimensionedScalar sigmaOm_;
	dimensionedScalar sigmaK_;
	dimensionedScalar ck_;
	dimensionedScalar Cles_;
	dimensionedScalar Cdmin_;
	dimensionedScalar Ct_;
	
	
	dimensionedScalar CL_;
	dimensionedScalar CEta_;
	dimensionedScalar Beta0_;
	
	
        // Fields
        volScalarField k_;
        volScalarField omega_;
	volScalarField nuSgs_;
	volSymmTensorField sgsstress;
	volScalarField epsilon_;
//	volScalarField epsilonH_;
	volScalarField alpha_;
//---------------------------------------------------------------------	
        volScalarField k_res_;
        volScalarField kM_av_;
        volScalarField k_tot_;
	volScalarField k_ratio_;
	volScalarField epsi_ratio_;
	volScalarField UU_av_;
        	
	volScalarField L_;
	volScalarField LI;
	volScalarField Leta;

	volScalarField Pn_;
	volScalarField Pt_;
	volScalarField Dn_;
        	
	volVectorField U_av_;
        	

	dimensionedScalar kr_;
         Switch restart_;

//----------------------------------------------------------------------
       // Switches for model control
	const Switch  pureRANS_;
        const Switch DurbinLimiter_;
	dimensionedScalar Ca_;
	dimensionedScalar f_PANS_;
	
	// unified stuff
	volScalarField transfer_;
        volScalarField Lrans_;
	dimensionedScalar lMin_;
	volScalarField filterW_;
	volScalarField tauL_;
	dimensionedScalar tauLMin_;
	dimensionedScalar epsiMin_;
	volScalarField Cmu_;
	volScalarField Ret;
	volScalarField ratio_;
	volScalarField tscale_;
	volScalarField Sk;
	volScalarField Ckd_;

// CES_KOS stuffs
        dimensionedScalar beta_eta;
        dimensionedScalar Tav_;
        volScalarField Ltot_;
        volScalarField cOmega2m_;
        volScalarField R_;
        volScalarField epsiM_av_;
        volScalarField epsi_res_;
        volScalarField epsi_tot_;
        volScalarField Lr_;
        volScalarField Lrimp_;
        volScalarField Lrimp2_;
        volScalarField Deltar_;
        volScalarField gradU2_av_;
	const Switch  L_act_;
	const Switch  L_imp_;
	const Switch  L_imp_new_;
	const Switch  R_PITM_;
	const Switch  R_PANS_;
	const Switch  R_PITM_actual_;
	autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Protected Member Functions


public:

    //- Runtime type information
    TypeName("CES_KOS");


    // Constructors

        //- Construct from components
        CES_KOS
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~CES_KOS()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }
        
        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                 new volScalarField("DkEff", nuSgs_/sigmaK_ + nu())
            );
        }

        //- Return the effective diffusivity for skL
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", nuSgs_/sigmaOm_ + nu())
            );
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }


        //- Return the turbulence kinetic energy dissipation rate      
         virtual tmp<volScalarField> epsilon() const
          {
             return tmp<volScalarField>
                (
                new volScalarField ("epsilon_", k_*omega_+sqr(Ct_)*nu()*sqr(omega_) )
                );
          }
        



        //- Return SGS viscosity
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const;

        //- Return the deviatoric part of the divergence of Beff
        //  i.e. the additional term in the filtered NSE.
        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Solve the turbulence equations (k-skL) and correct the turbulence
        //  viscosity
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HybridLesRasModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
