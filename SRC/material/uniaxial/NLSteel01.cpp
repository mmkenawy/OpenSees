/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.20 $
// $Date: 2008-08-26 16:35:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NLSteel01.cpp,v $
                                                                        
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// NLSteel01.
//
// What: "@(#) NLSteel01.C, revA"


#include <NLSteel01.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>

#include <math.h>
#include <float.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_NLSteel01()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[13];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial NLSteel01 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 9 && numData != 13) {
    opserr << "Invalid #args, want: uniaxialMaterial NLSteel01 " << iData[0] << " fy? E? b? <a1? a2? a3? a4?> epup? epun? ep0p? ep0n? maxpDam? maxnDam?" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial NLSteel01 " << iData[0] << " fy? E? b? <a1? a2? a3? a4?> epup? epun? ep0p? ep0n? maxpDam? maxnDam?" << endln;
    return 0;
  }

  if (numData == 9) {
    dData[9] = STEEL_01_DEFAULT_A1;
    dData[10] = STEEL_01_DEFAULT_A2;
    dData[11] = STEEL_01_DEFAULT_A3;
    dData[12] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new NLSteel01(iData[0], dData[0], dData[1],
			    dData[2], dData[3], dData[4], 
			    dData[5], dData[6], dData[7], dData[8],
				dData[9], dData[10], dData[11], dData[12]);

  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type NLSteel01 Material\n";
    return 0;
  }

  return theMaterial;
}



NLSteel01::NLSteel01
(int tag, double FY, double E, double B,
 double epup, double epun, double ep0p, double ep0n,
 double maxpDam, double maxnDam, double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_NLSteel01),
   fy(FY), E0(E), b(B), initpStrain(epup), initnStrain(epun),
   maxpStrain(ep0p), maxnStrain(ep0n),maxpDamage(maxpDam),
   maxnDamage(maxnDam), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
   // History variables
	Energy = 0;	//by SAJalali
	
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cnlstrain = 0.0;
   Custress = 0.0;
   Cstress = 0.0;
   Cpdamage = 0.0;
   Cndamage = 0.0;
   Cutangent = E0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tnlstrain = 0.0;
   Tustress = 0.0;
   Tstress = 0.0;
   Tpdamage = 0.0;
   Tndamage = 0.0;
   Tutangent = E0;
   Ttangent = E0;

//// AddingSensitivity:BEGIN /////////////////////////////////////
//	parameterID = 0;
//	SHVs = 0;
//// AddingSensitivity:END //////////////////////////////////////
}

NLSteel01::NLSteel01():UniaxialMaterial(0,MAT_TAG_NLSteel01),
 fy(0.0), E0(0.0), b(0.0), initpStrain(0.0), initnStrain(0.0),
 maxpStrain(0.0), maxnStrain(0.0), maxpDamage(0.0), maxnDamage(0.0),
 a1(0.0), a2(0.0), a3(0.0), a4(0.0)
{
	Energy = 0;	//by SAJalali

//// AddingSensitivity:BEGIN /////////////////////////////////////
//	parameterID = 0;
//	SHVs = 0;
//// AddingSensitivity:END //////////////////////////////////////

}

NLSteel01::~NLSteel01 ()
{
//// AddingSensitivity:BEGIN /////////////////////////////////////
//	if (SHVs != 0)
//		delete SHVs;
//// AddingSensitivity:END //////////////////////////////////////
}

int NLSteel01::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   //Tnlstrain = Cnlstrain; // this is updated in setNLStrain
   Tustress = Custress;
   Tstress = Cstress;
   Tpdamage = Cpdamage;
   Tndamage = Cndamage;
   Tutangent = Cutangent;
   Ttangent = Ctangent;

   // Update the damage parameters in both tension and compression
   //applyDamage();

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   return 0;
}

int NLSteel01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   //Tnlstrain = Cnlstrain; // this is updated in setNLStrain
   Tustress = Custress;
   Tstress = Cstress;
   Tpdamage = Cpdamage;
   Tndamage = Cndamage;
   Tutangent = Cutangent;
   Ttangent = Ctangent;

   // Update the damage parameters in both tension and compression
   //applyDamage();

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void NLSteel01::applyDamage (void)
{
  // hard-coded damage parameters for non-local degradation of the yield stress
  double m = 0.0; // (overly) non-local strain averaging parameter
//  double maxpDamage = 0.8; // maximum allowable + damage value < 1.0
//  double maxnDamage = 0.8; // maximum allowable - damage value < 1.0
//  double maxpStrain = 0.03; // + strain at which the material is fully damaged
//  double maxnStrain = 0.03; // - strain at which the material is fully damaged
//  double initpStrain = 0.01; // + strain at which damage initiates
//  double initnStrain = 0.01; // - strain at which damage initiates
    
  // pre-compute strain ranges over which damage is actively evolving
  double pRange = maxpStrain - initpStrain;
  double nRange = maxnStrain - initnStrain;

  // update the averaged (overly) non-local strain
  double nlstrain = m*Tnlstrain + (1.0-m)*Tstrain;

  // check which damage value needs to be updated depending on the non-local strain
  if (nlstrain > 0.0) {
    // update the damage parameter is tension
    Tpdamage = fmin(maxpDamage,fmax(Cpdamage,1.0-(maxpStrain-nlstrain)/pRange));
  } else {
    // update the damage parameter is compression
    Tndamage = fmin(maxnDamage,fmax(Cndamage,1.0-(maxnStrain+nlstrain)/nRange));
  }
}

void NLSteel01::determineTrialState (double dStrain)
{
      // compute the initial envelope stress intercept at zero strain;
      double fyOneMinusB = fy * (1.0 - b);

      // compute the hardening modulus: scale the initial elastic modulus
      // by the reduction factor b = Esh/E0 < 1.0
      double Esh = b*E0;

      // compute the strain at first yielding
      double epsy = fy/E0;

      // compute the unshifted stress on the envelope
      double c1 = Esh*Tstrain;

      // compute the positive envelope stress intercept at zero strain
      double c2 = TshiftN*fyOneMinusB;

      // compute the negative envelope stress intercept at zero strain
      double c3 = TshiftP*fyOneMinusB;

      // compute the trial elastic stress
      double c = Cstress + E0*dStrain;

      /**********************************************************
         removal of the following lines due to problems with
	 optimization may be required (e.g. on gnucc compiler
         with optimization turned on & -ffloat-store option not
         used) .. replace them with line that follows but which 
         now requires 2 function calls to achieve same result !!
      ************************************************************/
      // apply the damage
      //applyDamage();
      // compute the stress on the positive loading envelope,
      // scaled by the positive damage factor
      double c1c3 = (c1 + c3);
      // check for loading on the positive envelope
      if ((1.0 - Tpdamage)*c1c3 < c) {
	// Actively update the damage parameters only when loading on the envelopes
    applyDamage();
	// if the trial stress exceeds the positive envelope stress
	// project the stress back onto the positive envelope
	Tstress = (1.0 - Tpdamage)*c1c3;
	// the material is loading along the positive envelope,
	// use the hardening stiffness scaled by the positive damage factor
	Ttangent = (1.0 - Tpdamage)*Esh;
      } else {
	// otherwise, assume the trial stress is elastic
	// (will check for loading on negative envelope later on...)
	Tstress = c;
	// if the computed stress is the same as the trial (elastic) stress,
	// use the elastic stiffness
	Ttangent = E0;
      }

      // compute the stress on the negative loading envelope,
      // scaled by the negative damage factor
      double c1c2 = (c1 - c2);
      // check for loading on the negative envelope
      if ((1.0 - Tndamage)*c1c2 > Tstress) {
	// Actively update the damage parameters only when loading on the envelopes
	applyDamage();
	// if the trial stress is below the negative envelope stress
	// project the stress back onto the negative envelope
	Tstress = (1.0 - Tndamage)*c1c2;
	// the material is loading along the negative envelope,
	// use the hardening stiffness scaled by the negative damage factor
	Ttangent = (1.0 - Tndamage)*Esh;
      }

      /* ***********************************************************
      and replace them with:

      Tustress = fmax((c1-c2), fmin((c1+c3),c));
      **************************************************************/

      //
      // Determine if a load reversal has occurred due to the trial strain
      //

      // Determine initial loading condition
      if (Tloading == 0 && dStrain != 0.0) {
	  if (dStrain > 0.0)
	    Tloading = 1;
	  else
	    Tloading = -1;
      }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  // set the new maximum strain achieved by the material
	  if (Cstrain > TmaxStrain) {
	    TmaxStrain = Cstrain;
	  }
	  // grow the negative envelope due to hardening effects
	  TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  // set the new minimum strain achieved by the material
	  if (Cstrain < TminStrain) {
	    TminStrain = Cstrain;
	  }
	  // grow the positive envelope due to hardening effects
	  TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
      }
}

void NLSteel01::detectLoadReversal (double dStrain)
{
   // Determine initial loading condition
   if (Tloading == 0 && dStrain != 0.0)
   {
      if (dStrain > 0.0)
         Tloading = 1;
      else
         Tloading = -1;
   }

   double epsy = fy/E0;

   // Transition from loading to unloading, i.e. positive strain increment
   // to negative strain increment
   if (Tloading == 1 && dStrain < 0.0)
   {
      Tloading = -1;
      if (Cstrain > TmaxStrain)
         TmaxStrain = Cstrain;
      TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
   }

   // Transition from unloading to loading, i.e. negative strain increment
   // to positive strain increment
   if (Tloading == -1 && dStrain > 0.0)
   {
      Tloading = 1;
      if (Cstrain < TminStrain)
         TminStrain = Cstrain;
      TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
   }
}

int NLSteel01::setNLStrain(double nlstrain)
{
    Tnlstrain = nlstrain;
    return 0;
}

double NLSteel01::getStrain ()
{
   return Tstrain;
}

double NLSteel01::getStress ()
{
   return Tstress;
}

double NLSteel01::getTangent ()
{
   return Ttangent;
}

int NLSteel01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;

   // State variables
   //by SAJalali
   Energy += 0.5*(Tstress + Cstress)*(Tstrain - Cstrain);

   Cstrain = Tstrain;
   Cnlstrain = Tnlstrain;
   Cstress = Tstress;
   Custress = Tustress;
   Cpdamage = Tpdamage;
   Cndamage = Tndamage;
   Cutangent = Tutangent;
   Ctangent = Ttangent;

   return 0;
}

int NLSteel01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tnlstrain = Cnlstrain;
   Tstress = Cstress;
   Tustress = Custress;
   Tpdamage = Cpdamage;
   Tndamage = Cndamage;
   Tutangent = Cutangent;
   Ttangent = Ctangent;

   return 0;
}

int NLSteel01::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cnlstrain = 0.0;
   Cstress = 0.0;
   Custress = 0.0;
   Cpdamage = 0.0;
   Cndamage = 0.0;
   Cutangent = E0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tnlstrain = 0.0;
   Tstress = 0.0;
   Tustress = 0.0;
   Tpdamage = 0.0;
   Tndamage = 0.0;
   Tutangent = E0;
   Ttangent = E0;

//// AddingSensitivity:BEGIN /////////////////////////////////
//	if (SHVs != 0)
//		SHVs->Zero();
//// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* NLSteel01::getCopy ()
{
   NLSteel01* theCopy = new NLSteel01(this->getTag(), fy, E0, b,
				  initpStrain, initnStrain, maxpStrain, maxnStrain,
				  maxpDamage, maxnDamage, a1, a2, a3, a4);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cnlstrain = Cnlstrain;
   theCopy->Cstress = Cstress;
   theCopy->Custress = Custress;
   theCopy->Cpdamage = Cpdamage;
   theCopy->Cndamage = Cndamage;
   theCopy->Cutangent = Cutangent;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tnlstrain = Tnlstrain;
   theCopy->Tstress = Tstress;
   theCopy->Tustress = Tustress;
   theCopy->Tpdamage = Tpdamage;
   theCopy->Tndamage = Tndamage;
   theCopy->Tutangent = Tutangent;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int NLSteel01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(27);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = initpStrain;
   data(5) = initnStrain;
   data(6) = maxpStrain;
   data(7) = maxnStrain;
   data(8) = maxpDamage;
   data(9) = maxnDamage;
   data(10) = a1;
   data(11) = a2;
   data(12) = a3;
   data(13) = a4;

   // History variables from last converged state
   data(14) = CminStrain;
   data(15) = CmaxStrain;
   data(16) = CshiftP;
   data(17) = CshiftN;
   data(18) = Cloading;

   // State variables from last converged state
   data(19) = Cstrain;
   data(20) = Cnlstrain;
   data(21) = Cstress;
   data(22) = Custress;
   data(23) = Cpdamage;
   data(24) = Cndamage;
   data(25) = Cutangent;
   data(26) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "NLSteel01::sendSelf() - failed to send data\n";

   return res;
}

int NLSteel01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(27);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "NLSteel01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      initpStrain = data(4);
      initnStrain = data(5);
      maxpStrain = data(6);
      maxnStrain = data(7);
      maxpDamage = data(8);
      maxnDamage = data(9);
      a1 = data(10);
      a2 = data(11);
      a3 = data(12);
      a4 = data(13);

      // History variables from last converged state
      CminStrain = data(14);
      CmaxStrain = data(15);
      CshiftP = data(16);
      CshiftN = data(17);
      Cloading = int(data(18));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;

      // State variables from last converged state
      Cstrain = data(19);
      Cnlstrain = data(20);
      Cstress = data(21);
      Custress = data(22);
      Cpdamage = data(23);
      Cndamage = data(24);
      Cutangent = data(25);
      Ctangent = data(26);

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tnlstrain = Cnlstrain;
      Tstress = Cstress;
      Tustress = Custress;
      Tpdamage = Cpdamage;
      Tndamage = Cndamage;
      Tutangent = Cutangent;
      Ttangent = Ctangent;
   }
    
   return res;
}

void NLSteel01::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "NLSteel01 tag: " << this->getTag() << endln;
    s << "  fy: " << fy << " ";
    s << "  E0: " << E0 << " ";
    s << "   b: " << b << " ";
    s << "  a1: " << a1 << " ";
    s << "  a2: " << a2 << " ";
    s << "  a3: " << a3 << " ";
    s << "  a4: " << a4 << " ";
    s << "  initpStrain: " << initpStrain << " ";
    s << "  initnStrain: " << initnStrain << " ";
    s << "  maxpStrain: " << maxpStrain << " ";
    s << "  maxnStrain: " << maxnStrain << " ";
    s << "  maxpDamage: " << maxpDamage << " ";
    s << "  maxnDamage: " << maxnDamage << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"NLSteel01\", ";
	s << "\"E\": " << E0 << ", ";
	s << "\"fy\": " << fy << ", ";
    s << "\"b\": " << b << ", ";
    s << "\"a1\": " << a1 << ", ";
    s << "\"a2\": " << a2 << ", ";
    s << "\"a3\": " << a3 << ", ";
    s << "\"a4\": " << a4 << ", ";
    s << "\"initpStrain\": " << initpStrain << ", ";
    s << "\"initnStrain\": " << initnStrain << ", ";
    s << "\"maxpStrain\": " << maxpStrain << ", ";
    s << "\"maxnStrain\": " << maxnStrain << ", ";
    s << "\"maxpDamage\": " << maxpDamage << ", ";
    s << "\"maxnDamage\": " << maxnDamage << "}";
  }
  
}




//// AddingSensitivity:BEGIN ///////////////////////////////////
//int
//NLSteel01::setParameter(const char **argv, int argc, Parameter &param)
//{
//
//  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
//    param.setValue(fy);
//    return param.addObject(1, this);
//  }
//  if (strcmp(argv[0],"E") == 0) {
//    param.setValue(E0);
//    return param.addObject(2, this);
//  }
//  if (strcmp(argv[0],"b") == 0) {
//    param.setValue(b);
//    return param.addObject(3, this);
//  }
//  if (strcmp(argv[0],"a1") == 0) {
//    param.setValue(a1);
//    return param.addObject(4, this);
//  }
//  if (strcmp(argv[0],"a2") == 0) {
//    param.setValue(a2);
//    return param.addObject(5, this);
//  }
//  if (strcmp(argv[0],"a3") == 0) {
//    param.setValue(a3);
//    return param.addObject(6, this);
//  }
//  if (strcmp(argv[0],"a4") == 0) {
//    param.setValue(a4);
//    return param.addObject(7, this);
//  }
//
//  return -1;
//}
//
//
//
//int
//NLSteel01::updateParameter(int parameterID, Information &info)
//{
//	switch (parameterID) {
//	case -1:
//		return -1;
//	case 1:
//		this->fy = info.theDouble;
//		break;
//	case 2:
//		this->E0 = info.theDouble;
//		break;
//	case 3:
//		this->b = info.theDouble;
//		break;
//	case 4:
//		this->a1 = info.theDouble;
//		break;
//	case 5:
//		this->a2 = info.theDouble;
//		break;
//	case 6:
//		this->a3 = info.theDouble;
//		break;
//	case 7:
//		this->a4 = info.theDouble;
//		break;
//	default:
//		return -1;
//	}
//
//	Tutangent = E0;          // Initial stiffness
//
//	return 0;
//}
//
//
//
//
//int
//NLSteel01::activateParameter(int passedParameterID)
//{
//	parameterID = passedParameterID;
//
//	return 0;
//}
//
//
//
//double
//NLSteel01::getStressSensitivity(int gradIndex, bool conditional)
//{
//	// Initialize return value
//	double gradient = 0.0;
//
//
//	// Pick up sensitivity history variables
//	double CstrainSensitivity = 0.0;
//	double CstressSensitivity = 0.0;
//	if (SHVs != 0) {
//		CstrainSensitivity = (*SHVs)(0,gradIndex);
//		CstressSensitivity = (*SHVs)(1,gradIndex);
//	}
//
//
//	// Assign values to parameter derivatives (depending on what's random)
//	double fySensitivity = 0.0;
//	double E0Sensitivity = 0.0;
//	double bSensitivity = 0.0;
//	if (parameterID == 1) {
//		fySensitivity = 1.0;
//	}
//	else if (parameterID == 2) {
//		E0Sensitivity = 1.0;
//	}
//	else if (parameterID == 3) {
//		bSensitivity = 1.0;
//	}
//
//
//	// Compute min and max stress
//	double Tstress;
//	double dStrain = Tstrain-Cstrain;
//	double sigmaElastic = Cstress + E0*dStrain;
//	double fyOneMinusB = fy * (1.0 - b);
//	double Esh = b*E0;
//	double c1 = Esh*Tstrain;
//	double c2 = TshiftN*fyOneMinusB;
//	double c3 = TshiftP*fyOneMinusB;
//	double sigmaMax = c1+c3;
//	double sigmaMin = c1-c2;
//
//
//	// Evaluate stress sensitivity
//	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
//		Tstress = sigmaMax;
//		gradient = E0Sensitivity*b*Tstrain
//				 + E0*bSensitivity*Tstrain
//				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
//	}
//	else {
//		Tstress = sigmaElastic;
//		gradient = CstressSensitivity
//			     + E0Sensitivity*(Tstrain-Cstrain)
//				 - E0*CstrainSensitivity;
//	}
//	if (sigmaMin > Tstress) {
//		gradient = E0Sensitivity*b*Tstrain
//			     + E0*bSensitivity*Tstrain
//				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
//	}
//
//	return gradient;
//}
//
//
//
//
//double
//NLSteel01::getInitialTangentSensitivity(int gradIndex)
//{
//	// For now, assume that this is only called for initial stiffness
//	if (parameterID == 2) {
//		return 1.0;
//	}
//	else {
//		return 0.0;
//	}
//}
//
//
//int
//NLSteel01::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
//{
//	if (SHVs == 0) {
//		SHVs = new Matrix(2,numGrads);
//	}
//
//
//	// Initialize unconditaional stress sensitivity
//	double gradient = 0.0;
//
//
//	// Pick up sensitivity history variables
//	double CstrainSensitivity = 0.0;
//	double CstressSensitivity	 = 0.0;
//	if (SHVs != 0) {
//		CstrainSensitivity = (*SHVs)(0,gradIndex);
//		CstressSensitivity = (*SHVs)(1,gradIndex);
//	}
//
//
//	// Assign values to parameter derivatives (depending on what's random)
//	double fySensitivity = 0.0;
//	double E0Sensitivity = 0.0;
//	double bSensitivity = 0.0;
//	if (parameterID == 1) {
//		fySensitivity = 1.0;
//	}
//	else if (parameterID == 2) {
//		E0Sensitivity = 1.0;
//	}
//	else if (parameterID == 3) {
//		bSensitivity = 1.0;
//	}
//
//
//	// Compute min and max stress
//	double Tstress;
//	double dStrain = Tstrain-Cstrain;
//	double sigmaElastic = Cstress + E0*dStrain;
//	double fyOneMinusB = fy * (1.0 - b);
//	double Esh = b*E0;
//	double c1 = Esh*Tstrain;
//	double c2 = TshiftN*fyOneMinusB;
//	double c3 = TshiftP*fyOneMinusB;
//	double sigmaMax = c1+c3;
//	double sigmaMin = c1-c2;
//
//
//	// Evaluate stress sensitivity ('gradient')
//	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
//		Tstress = sigmaMax;
//		gradient = E0Sensitivity*b*Tstrain
//				 + E0*bSensitivity*Tstrain
//				 + E0*b*TstrainSensitivity
//				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
//	}
//	else {
//		Tstress = sigmaElastic;
//		gradient = CstressSensitivity
//			     + E0Sensitivity*(Tstrain-Cstrain)
//				 + E0*(TstrainSensitivity-CstrainSensitivity);
//	}
//	if (sigmaMin > Tstress) {
//		gradient = E0Sensitivity*b*Tstrain
//			     + E0*bSensitivity*Tstrain
//			     + E0*b*TstrainSensitivity
//				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
//	}
//
//
//	// Commit history variables
//	(*SHVs)(0,gradIndex) = TstrainSensitivity;
//	(*SHVs)(1,gradIndex) = gradient;
//
//	return 0;
//}
//
//// AddingSensitivity:END /////////////////////////////////////////////

