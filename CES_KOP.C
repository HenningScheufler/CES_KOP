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

#include "CES_KOS.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace HybridLesRasModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(CES_KOS, 0);
addToRunTimeSelectionTable(HybridLesRasModel, CES_KOS, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CES_KOS::CES_KOS
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    HybridLesRasModel(modelName, U, phi, transport, turbulenceModelName),

    cOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cOmega1",
            coeffDict_,
            0.49
        )
    ),
    cOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cOmega2",
            coeffDict_,
            0.072
        )
    ),
    cOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cOmega",
            coeffDict_,
            1.1
        )
    ),
    sigmaOm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOm",
            coeffDict_,
            1.8
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),
    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.09
        )
    ),
    Cles_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cles",
            coeffDict_,
            1.0
        )
    ),

    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.16
        )
    ),
    CEta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            11.0
        )
    ),

     Ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct",
            coeffDict_,
            0.3
        )
    ),


    Cdmin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdmin",
            coeffDict_,

           -0.05
        )
    ),
    pureRANS_(coeffDict_.lookupOrAddDefault<Switch>("pureRANS", false)),
    L_act_(coeffDict_.lookupOrAddDefault<Switch>("L_act", true)),
    L_imp_(coeffDict_.lookupOrAddDefault<Switch>("L_imp", true)),
    L_imp_new_(coeffDict_.lookupOrAddDefault<Switch>("L_imp_new", true)),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

     epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
     dimensionedScalar("epsilon",dimensionSet(0, 2, -3, 0, 0, 0, 0),SMALL)
    ),


   alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
     dimensionedScalar("alpha",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),

// Moving time averages

 U_av_
    (
        IOobject
        (
            "U_av",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("U_av",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    ),
 UU_av_
    (
        IOobject
        (
            "UU_av",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("UU_av",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    ),

    k_res_
    (
        IOobject
        (
            "k_res",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("k_res",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    ),

   k_tot_
    (
        IOobject
        (
            "k_tot",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("k_tot",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    ),

    k_ratio_
    (
        IOobject
        (
            "k_ratio",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("k_ratio",dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
    ),
epsi_ratio_
    (
        IOobject
        (
            "epsi_ratio",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsi_ratio",dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
    ),


    kM_av_
    (
        IOobject
        (
            "kM_av",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kM_av",dimensionSet(0, 2, -2, 0, 0, 0, 0),0.0)
    ),
    kr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kr",
            coeffDict_,
            0.1
        )
    ),



    restart_(coeffDict_.lookupOrAddDefault<Switch>("restart", true)),
 
   nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
    transfer_
    (
        IOobject
        (
            "transfer",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("transfer",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),

  cOmega2m_
    (
        IOobject
        (
            "cOmega2m",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cOmega2m",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),

  R_
    (
        IOobject
        (
            "R",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("R",dimensionSet(0, 0, 0, 0, 0, 0, 0),1.0)
    ),


    Lrans_
    (
        IOobject
        (
            "Lrans",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Lrans",dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL)
    ),
    lMin_("lMin", dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL),
    tauLMin_("tauLMin", dimensionSet(0, 0, 1, 0, 0, 0, 0),SMALL),
    epsiMin_("epsiMin", dimensionSet(0, 2, -3, 0, 0, 0, 0),SMALL),

    
    filterW_
    (
        IOobject
        (
            "filterW",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("filterW",dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL)
    ),
    tauL_
    (
        IOobject
        (
            "tauL",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("tauL",dimensionSet(0, 0, 1, 0, 0, 0, 0),SMALL)
    ), 
    tscale_
    (
        IOobject
        (
            "tscale",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("tscale",dimensionSet(0, 0, 1, 0, 0, 0, 0),SMALL)
    ), 
    Cmu_
    (
        IOobject
        (
            "Cmu",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Cmu",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),
    Ret
    (
        IOobject
        (
            "Ret",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Ret",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),
    Sk
    (
        IOobject
        (
            "Sk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sk",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),
    ratio_
    (
        IOobject
        (
            "ratio",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("ratio",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),
    sgsstress
    (
      IOobject
      (
        "sgsstress",
        runTime_.timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedSymmTensor("sgsstress",dimensionSet(0, 2, -2, 0, 0, 0, 0),symmTensor::zero)
    ),
    Ckd_
    (
       IOobject
       (
          "Ckd",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("Ckd",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),

       Pn_
    (
       IOobject
       (
          "Pn",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("Pn",dimensionSet(0, 2, -3, 0, 0, 0, 0),SMALL)
    ),

  Dn_
    (
       IOobject
       (
          "Dn",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("Dn",dimensionSet(0, 2, -3, 0, 0, 0, 0),SMALL)
    ),


  Pt_
    (
       IOobject
       (
          "Pt",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("Pt",dimensionSet(0, 2, -3, 0, 0, 0, 0),SMALL)
    ),



    L_
    (
       IOobject
       (
          "L",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("L",dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL)
    ),

    LI
    (
       IOobject
       (
          "LI",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("LI",dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL)
    ),

   Leta
    (
       IOobject
       (
          "Leta",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("Leta",dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL)
    ),

     f_PANS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "f_PANS",
            coeffDict_,
            0.4
        )
    ),

DurbinLimiter_(coeffDict_.lookupOrAddDefault<Switch>("DurbinLimiter", false)),
R_PITM_(coeffDict_.lookupOrAddDefault<Switch>("R_PITM", false)),
R_PANS_(coeffDict_.lookupOrAddDefault<Switch>("R_PANS`", false)),
R_PITM_actual_(coeffDict_.lookupOrAddDefault<Switch>("R_PITM_actual", false)),

  Ca_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ca",
            coeffDict_,
            0.4
        )
    ),
Beta0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Beta0",
            coeffDict_,
            0.44
        )
    ),


Tav_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Tav",
            coeffDict_,
            0.0
        )
    ),

 beta_eta
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta_eta",
            coeffDict_,
            0.026
        )
    ),
   Ltot_
    (
       IOobject
       (
          "Ltot",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("Ltot",dimensionSet(0, 1, 0, 0, 0, 0, 0),SMALL)
    ),

    Lr_
    (
        IOobject
        (
            "Lr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Lr",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),
Lrimp_
    (
        IOobject
        (
            "Lrimp",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Lrimp",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),

Lrimp2_
    (
        IOobject
        (
            "Lrimp2",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Lrimp2",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),

 Deltar_
    (
        IOobject
        (
            "Deltar",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Deltar",dimensionSet(0, 0, 0, 0, 0, 0, 0),SMALL)
    ),


 epsi_res_
    (
        IOobject
        (
            "epsi_res",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsi_res",dimensionSet(0, 2, -3, 0, 0, 0, 0),0.0)
    ),
    epsiM_av_
    (
        IOobject
        (
            "epsiM_av",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsiM_av",dimensionSet(0, 2, -3, 0, 0, 0, 0),0.0)
    ),

   epsi_tot_ 
    (
        IOobject
        (
            "epsi_tot",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsi_tot",dimensionSet(0, 2, -3, 0, 0, 0, 0),0.0)
    ),
gradU2_av_
    (
        IOobject
        (
            "gradU2_av",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gradU2_av",dimensionSet(0, 0, -2, 0, 0, 0, 0),0.0)
    ),


    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);
    Info << "kMin_=" << kMin_ << "omegaMin=" << omegaMin_ << endl;
    
    Info << ">>>>Tav=    " <<Tav_ <<endl;
    Info << ">>>>kr=     " <<kr_ <<endl;
    Info << ">>>>beta_eta=     " <<beta_eta <<endl;
    // unified time scale calculation
    filterW_ = delta();
    Lrans_ = sqrt(k_)/(omega_ + omegaMin_) + lMin_;

// initialize the Le with Lrans
    Ltot_= Lrans_;

    tscale_ = 1./(omega_ + omegaMin_)+ tauLMin_;
    ratio_ = (filterW_ / sqrt(k_ + kMin_) / Cles_) / tscale_;

    tauL_ = tscale_ ;
    
    Cmu_=0.09 * pow(alpha_, 3.0); 
    Info << ">>>>>>>>>>>>>>>>>>>>>>Cmu above >>>>>>>>>>>>>>>" <<endl;
    nuSgs_ = Cmu_ * k_ * tauL_;
    nuSgs_.correctBoundaryConditions();
    
    // Turbulence Reynolds number
    Ret = 
    (
      nuSgs_ / nu()
    );
    

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CES_KOS::correct(const tmp<volTensorField>& gradU)
{
    HybridLesRasModel::correct(gradU);
    
    // S_ii and Omega_ij
    volSymmTensorField Sijd=dev(symm(gradU()));
    volTensorField Rotij=skew(-gradU());
    
    // S2
    volScalarField S2 = magSqr(Sijd);

    // generation term
    volScalarField G
    (
        "HybridLesRasModel::G",
        nuSgs_*2.*S2
    );
//
// C_av= 1.0/ (1 + (dt/(Tav-dt)));
 scalar C_av_ = 1./(1. + runTime_.deltaT().value()/(Tav_.value()-runTime_.deltaT().value()));
    Info << ">>>>C_av=     " <<C_av_ <<endl;
 dimensionedScalar Beta=cOmega2_/(ck_*cOmega1_);
    Info<<">>>Beta=  "<< Beta.value()<<endl;

  // Advance moving time averages
    if (restart_)
    {
      U_av_        = U();
      UU_av_       = tr( (U() * U()) );
      kM_av_       = k_;
      epsiM_av_    = epsilon();
      gradU2_av_   = gradU() && gradU();
      restart_ = false;
      Info << restart_<<" >>>>Restarted the moving time averages" << endl;
    }
  else
    {
      U_av_        = (1. - C_av_) * U() + C_av_ * U_av_;
      UU_av_       = (1. - C_av_) * tr( (U_ * U_) ) + C_av_ * UU_av_;
      kM_av_       = (1. - C_av_) * k_ + C_av_ * kM_av_;
      epsiM_av_      = (1. - C_av_) * epsilon() + C_av_ * epsiM_av_;
      gradU2_av_ = (1. - C_av_) * (gradU() && gradU()) + C_av_ * gradU2_av_;

      Info << "restart_"<<restart_ << endl;
    }
    // Calculate relevant variables
     k_res_   = 0.5*(UU_av_ - tr( (U_av_ * U_av_) ));
     k_tot_   = kM_av_ + k_res_;
    epsi_res_ = nu() * (gradU2_av_ - ((fvc::grad(U_av_,"grad(U)")) && (fvc::grad(U_av_,"grad(U)"))));
    epsi_res_.max(1.e-15);

     epsi_tot_   = epsiM_av_ + epsi_res_;
     k_ratio_ = kM_av_ / ( k_tot_+ kMin_ );
     epsi_ratio_ = epsiM_av_ / ( epsi_tot_+ epsiMin_ );

    
    Ltot_ = pow(k_tot_,1.5)/(epsi_tot_+ epsiMin_);
    //Lrans_ = sqrt(k_+kMin_)/(omega_ + omegaMin_) + lMin_;
    Lrans_ = pow(kM_av_,1.5)/(epsiM_av_ + epsiMin_) ;
    
    Lr_=Lrans_/(Ltot_+lMin_);
    Deltar_= filterW_/(Ltot_+lMin_);


    Lrimp2_= pow( 1 + pow(Deltar_,-3.0) , -1.0/3.0);


    forAll(Lrimp_,celli)
      {
      if(Deltar_[celli] < 1.0)
       {
        Lrimp_[celli]= (sqrt( 1 + 4*Beta.value()*(Beta.value()-1.0)*sqr(Deltar_[celli]))- 1.0) / ( 2.0* (Beta.value()-1.0)*(Deltar_[celli]) );
       }
      else
       {
       Lrimp_[celli]=1.0;
       }
      }
     


    if (L_act_)  {  
                   R_ = sqr(Lr_);
                   Info<< "L_act_ = true " <<endl;
                 }
  
    else if (L_imp_)
               {

                  if (L_imp_new_)
                     {
                      R_ = sqr(Lrimp2_);
                      Info<< "L_imp_ = true, L_imp_new " <<endl;

                     }
                  else
                    {
                      R_ = sqr(Lrimp_);
                      Info<< "L_imp_ = true, L_imp_old " <<endl;
                    }
                }

 
   if (R_PITM_) {
           
                   Info<< "R_PITM_ = true " <<endl;
                   if  (R_PITM_actual_){
                                         Info<< "R_PITM_actuual = true " <<endl;
                                         R_ = k_ratio_;
                                       }
                   else{  
                         Info<< "R_PITM_model = true " <<endl;
                         R_= pow(Lrimp2_,2.0/3.0);
                        }
                 }

    if (R_PANS_) {
                   Info<< "R_PANS_ = true " <<endl;
                   R_ = f_PANS_;
                 }



                cOmega2m_ = cOmega1_ + R_ * (cOmega2_/ck_ - cOmega1_) ;
     

    // Update omega at the wall
    #include "include/omegaWall.H"
    
    // omega-equation Cross diffusion term
    volScalarField cdomega=(fvc::grad(k_) & fvc::grad(omega_))/(k_+kMin_);

   // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
        cOmega1_*omega_*G/(k_+ kMin_)
      - fvm::Sp(cOmega2m_* omega_, omega_)    // dividing by ck comes from the transformation omega = ck * omega^* where omega^* is Bredbergs omega
      + fvm::SuSp
        (
            cOmega_*(nu()+nuSgs_)*cdomega/(omega_+omegaMin_),
            omega_
        )
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
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G       
      - fvm::Sp(1./tauL_, k_)  
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
    
    // Calculate test filtered fields for dynamic model coeffcient
    volSymmTensorField Lijd = (filter_(sqr(U())) - (sqr(filter_(U()))));
    volScalarField ktest = 0.5*tr(Lijd);
    Lijd = dev(Lijd);
    
    // Recalculate time scale
    tscale_ = 1./(omega_ + omegaMin_)+ tauLMin_;
//    ratio_ = (filterW_ / sqrt(k_ + kMin_) / Cles_) / tscale_;
//    transfer_ = min(ratio_,1.0);
    
    tauL_ = tscale_ ;
        
    



      epsilon_=k_*omega_ + sqr(Ct_)*nu()*sqr(omega_);

     LI = CL_* pow(k_,1.5)/epsilon_;

     Leta = CEta_ * pow(pow(nu(),3.0)/epsilon_,0.25);

     L_= LI;

    forAll(L_,celli)
    {
      if (k_ratio_[celli] > kr_.value()) // 
      {
     L_[celli]=max(LI[celli], Leta[celli]); 

      }
    }

     tmp<fvScalarMatrix> alphaEqn
     (
         fvm::laplacian(alpha_)
      ==
         fvm::Sp(1.0/sqr(L_), alpha_)
         - scalar(1.0)/( sqr(L_) )
     );
     alphaEqn().relax();
     solve(alphaEqn);
    // bound alpha
     alpha_.min(1.0);
     alpha_.max(0.0);

    Info << ":::::::::::Cmu down:::::" <<endl;
     Cmu_= 0.09 * pow(alpha_, 3.0);


    // Use Damping only in RANS region: use dynamic coeffcient in LES region
    nuSgs_ = Cmu_ * k_ * tauL_;
    
    nuSgs_.correctBoundaryConditions();
    
    // turbulence Reynolds number
    Ret = 
    (
      nuSgs_ / nu()
    );
    
    sgsstress = B();


}

tmp<volSymmTensorField> CES_KOS::B() const
{
    return ((2.0/3.0)*I)*k() - nuSgs()*twoSymm(fvc::grad(U()));
}


tmp<volSymmTensorField> CES_KOS::devBeff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> CES_KOS::divDevBeff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U) 
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


bool CES_KOS::read()
{
    if (HybridLesRasModel::read())
    {
        cOmega1_.readIfPresent(coeffDict());
        Beta0_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        CEta_.readIfPresent(coeffDict());
        cOmega2_.readIfPresent(coeffDict());
        cOmega_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaOm_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());
        Cles_.readIfPresent(coeffDict());
	Cdmin_.readIfPresent(coeffDict());
	Ca_.readIfPresent(coeffDict());
	Tav_.readIfPresent(coeffDict());
	beta_eta.readIfPresent(coeffDict());
	kr_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HybridLesRasModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
