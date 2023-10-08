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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel01.cpp,v $
                                                                        
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// Steel01. 
//
// What: "@(#) Steel01.C, revA"


#include <Steel01.h>
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
#include <iostream>


void *
OPS_Steel01()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01 " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01 " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    dData[3] = STEEL_01_DEFAULT_A1;
    dData[4] = STEEL_01_DEFAULT_A2;
    dData[5] = STEEL_01_DEFAULT_A3;
    dData[6] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new Steel01(iData[0], dData[0], dData[1], 
			    dData[2], dData[3], dData[4], 
			    dData[5], dData[6]);

  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel01 Material\n";
    return 0;
  }

  return theMaterial;
}



Steel01::Steel01
(int tag, double FY, double E, double B,
 double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_Steel01),
   fy(FY), E0(E), b(B), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
   // History variables
	Energy = 0;	//by SAJalali

   // State variables
   Cstrain   = 0.0; // strain
   Cnlstrain = 0.0; // non-local strain
   Cstress   = 0.0; // stress
   Ctangent  =  E0; // material stiffness (tangent)
   Cpstrain  = 0.0; // total plastic strain
   Ceps      = 0.0; // equivalent plastic strain
   Cdamage   = 0.0; // damage variable controlling compressive softening     

   Tstrain   = 0.0; // strain
   Tnlstrain = 0.0; // non-local strain
   Tstress   = 0.0; // stress
   Ttangent  =  E0; // material stiffness (tangent)
   Tpstrain  = 0.0; // total plastic strain
   Teps      = 0.0; // equivalent plastic strain
   Tdamage   = 0.0; // damage variable controlling compressive softening

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Steel01::Steel01():UniaxialMaterial(0,MAT_TAG_Steel01),
 fy(0.0), E0(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
{
	Energy = 0;	//by SAJalali

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

}

Steel01::~Steel01 ()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int Steel01::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
   Tstrain   = Cstrain;   // strain
   //Tnlstrain = Cnlstrain; // non-local strain (this is updated in setNLStrain)
   Tstress   = Cstress;   // stress
   Ttangent  = Ctangent;  // material stiffness (tangent)
   Tpstrain  = Cpstrain;  // total plastic strain
   Teps      = Ceps;      // equivalent plastic strain
   Tdamage   = Cdamage;   // damage variable controlling compressive softening

   // Set trial strain
   Tstrain = strain;

   // Calculate the trial state given the trial strain
   determineTrialState ();

   return 0;
}

int Steel01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   Tstrain   = Cstrain;   // strain
   //Tnlstrain = Cnlstrain; // non-local strain (this is updated in setNLStrain)
   Tstress   = Cstress;   // stress
   Ttangent  = Ctangent;  // material stiffness (tangent)
   Tpstrain  = Cpstrain;  // total plastic strain
   Teps      = Ceps;      // equivalent plastic strain
   Tdamage   = Cdamage;   // damage variable controlling compressive softening

   // Set trial strain
   Tstrain = strain;

   // Calculate the trial state given the trial strain
   determineTrialState ();

   stress  = Tstress;
   tangent = Ttangent;

   return 0;
}

//void Steel01::applyDamage (void)
//{
//  // hard-coded damage parameters for non-local degradation of the yield stress
//  double m = 1.5; // (overly) non-local strain averaging parameter
//  double maxpDamage = 0.8; // maximum allowable + damage value < 1.0
//  double maxnDamage = 0.8; // maximum allowable - damage value < 1.0
//  double maxpStrain = 0.03; // + strain at which the material is fully damaged
//  double maxnStrain = 0.03; // - strain at which the material is fully damaged
//  double initpStrain = 0.01; // + strain at which damage initiates
//  double initnStrain = 0.01; // - strain at which damage initiates
//    
//  // pre-compute strain ranges over which damage is actively evolving
//  double pRange = maxpStrain - initpStrain;
//  double nRange = maxnStrain - initnStrain;
//
//  // update the averaged (overly) non-local strain
//  double nlstrain = m*Tnlstrain + (1.0-m)*Tstrain;
//
//  // check which damage value needs to be updated depending on the non-local strain
//  if (nlstrain > 0.0) {
//    // update the damage parameter is tension
//    Tpdamage = fmin(maxpDamage,fmax(Cpdamage,1.0-(maxpStrain-nlstrain)/pRange));
//  } else {
//    // update the damage parameter is compression
//    Tndamage = fmin(maxnDamage,fmax(Cndamage,1.0-(maxnStrain+nlstrain)/nRange));
//  }
//}

void Steel01::determineTrialState ()
{
  // hard-coded damage parameters
  double m          = 1.0;  // (overly) non-local strain averaging parameter
  double initStrain = 0.0001; // (compressive) strain at which damage initiates
  double maxStrain  = 0.07; // 0.07 (compressive) strain at which material is fully damaged
  double maxDamage  = 1.0;  // maximum (compressive) damage
  
  //opserr << "Running Brian's model" << endln;

  // update the averaged (overly) non-local strain
  double nlstrain = m*Tnlstrain + (1.0-m)*Tstrain;
  
  // update the damage based on the non-local strain
  //Tdamage = std::min(std::max(Cdamage,tanh(-(nlstrain + initStrain)/(maxStrain - initStrain))),maxDamage);
  Tdamage = std::min(std::max(Cdamage,tanh(-(Tstrain + initStrain)/(maxStrain - initStrain))),maxDamage);
  //Tdamage = std::min(std::max(Cdamage,-(nlstrain + initStrain)/(maxStrain - initStrain)),maxDamage);

  // compute undamaged trial stress
  Tstress = E0*(Tstrain - Cpstrain);

  // check for violation of the yield constraint
  double phi = abs(Tstress - b*E0*Cpstrain) - fy;
  if (phi > 0.0) {
	// update the yield strain:

	// compute the loading direction
	double ny = copysign(1.0,Tstress - b*E0*Cpstrain);

	// compute the plastic strain increment
	double dep = phi/(E0*ny+b*E0*ny);

	// update the total and equivalent plastic strains
	Tpstrain = Cpstrain +     dep;
	Teps     = Ceps     + abs(dep);

	// update the undamaged stress
	Tstress = E0*(Tstrain - Tpstrain);
  }

  // Apply the damage scale factor to the stress
  if (Tstress < 0.0) Tstress *= (1.0 - Tdamage);

  /*// compute the non-local buckling strain
  //double bstrain = Tdamage*std::min(nlstrain - Cpstrain,0.0);
  double bstrain = Tdamage*std::min(Tstrain - Cpstrain,0.0);
  double fac = 1.0+0.8*std::min(copysign(Tdamage,Tstrain - Cpstrain),0.0);

  // compute the trial stress
  Tstress = E0*(Tstrain - Cpstrain - bstrain);

  // check for violation of the yield constrain
  double phi = abs(Tstress) - (fy + b*E0*Ceps);
  if (phi > 0.0) {
    // update the yield strain:

    // compute the loading direction
    double ny = copysign(1.0,Tstress);

    // compute the plastic strain increment
    double dep;
    //if (Tnlstrain < Cpstrain) {
    if (Tstrain < Cpstrain) {
      dep = phi/((1.0-Tdamage)*E0*ny+b*E0);
    } else {
      dep = phi/(E0*ny+b*E0);
    }

    // update the total and equivalent plastic strains
    Tpstrain = Cpstrain +     dep;
    Teps     = Ceps     + abs(dep);

    // update the buckling strain
    //bstrain = Tdamage*std::min(nlstrain - Tpstrain,0.0);
    bstrain = Tdamage*std::min(Tstrain - Tpstrain,0.0);

    // update the stress
    Tstress = E0*(Tstrain - Tpstrain - bstrain);
  }*/
  
}

int Steel01::setNLStrain(double nlstrain)
{
    Tnlstrain = nlstrain;
    return 0;
}

double Steel01::getStrain ()
{
   return Tstrain;
}

double Steel01::getStress ()
{
   return Tstress;
}

double Steel01::getTangent ()
{
   return Ttangent;
}

int Steel01::commitState ()
{
   // State variables
   //by SAJalali
   Energy += 0.5*(Tstress + Cstress)*(Tstrain - Cstrain);

   Cstrain   = Tstrain;   // strain
   Cnlstrain = Tnlstrain; // non-local strain
   Cstress   = Tstress;   // stress
   Ctangent  = Ttangent;  // material stiffness (tangent)
   Cpstrain  = Tpstrain;  // total plastic strain
   Ceps      = Teps;      // equivalent plastic strain
   Cdamage   = Tdamage;   // damage variable controlling compressive softening

   return 0;
}

int Steel01::revertToLastCommit ()
{
   // Reset trial state variables to last committed state
   Tstrain   = Cstrain;   // strain
   Tnlstrain = Cnlstrain; // non-local strain
   Tstress   = Cstress;   // stress
   Ttangent  = Ctangent;  // material stiffness (tangent)
   Tpstrain  = Cpstrain;  // total plastic strain
   Teps      = Ceps;      // equivalent plastic strain
   Tdamage   = Cdamage;   // damage variable controlling compressive softening

   return 0;
}

int Steel01::revertToStart ()
{
   // State variables
   Cstrain   = 0.0; // strain
   Cnlstrain = 0.0; // non-local strain
   Cstress   = 0.0; // stress
   Ctangent  =  E0; // material stiffness (tangent)
   Cpstrain  = 0.0; // total plastic strain
   Ceps      = 0.0; // equivalent plastic strain
   Cdamage   = 0.0; // damage variable controlling compressive softening     

   Tstrain   = 0.0; // strain
   Tnlstrain = 0.0; // non-local strain
   Tstress   = 0.0; // stress
   Ttangent  =  E0; // material stiffness (tangent)
   Tpstrain  = 0.0; // total plastic strain
   Teps      = 0.0; // equivalent plastic strain
   Tdamage   = 0.0; // damage variable controlling compressive softening

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* Steel01::getCopy ()
{
   Steel01* theCopy = new Steel01(this->getTag(), fy, E0, b,
				  a1, a2, a3, a4);

   // Converged state variables
   theCopy->Cstrain   = Cstrain;   // strain
   theCopy->Cnlstrain = Cnlstrain; // non-local strain
   theCopy->Cstress   = Cstress;   // stress
   theCopy->Ctangent  = Ctangent;  // material stiffness (tangent)
   theCopy->Cpstrain  = Cpstrain;  // total plastic strain
   theCopy->Ceps      = Ceps;      // equivalent plastic strain
   theCopy->Cdamage   = Cdamage;   // damage variable controlling compressive softening

   // Trial state variables
   theCopy->Tstrain   = Tstrain;   // strain
   theCopy->Tnlstrain = Tnlstrain; // non-local strain
   theCopy->Tstress   = Tstress;   // stress
   theCopy->Ttangent  = Ttangent;  // material stiffness (tangent)
   theCopy->Tpstrain  = Tpstrain;  // total plastic strain
   theCopy->Teps      = Teps;      // equivalent plastic strain
   theCopy->Tdamage   = Tdamage;   // damage variable controlling compressive softening

   return theCopy;
}

int Steel01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(14);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = a1;
   data(5) = a2;
   data(6) = a3;
   data(7) = a4;

   // State variables from last converged state
   data( 8) = Cstrain  ; // strain
   data( 9) = Cnlstrain; // non-local strain
   data(10) = Cstress  ; // stress
   data(11) = Ctangent ; // material stiffness (tangent)
   data(12) = Cpstrain ; // total plastic strain
   data(13) = Ceps     ; // equivalent plastic strain
   data(14) = Cdamage  ; // damage variable controlling compressive softening  

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Steel01::sendSelf() - failed to send data\n";

   return res;
}

int Steel01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(14);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "Steel01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b  = data(3);
      a1 = data(4);
      a2 = data(5);
      a3 = data(6);
      a4 = data(7);

      // State variables from last converged state
      Cstrain   = data( 8); // strain
      Cnlstrain = data( 9); // non-local strain
      Cstress   = data(10); // stress
      Ctangent  = data(11); // material stiffness (tangent)
      Cpstrain  = data(12); // total plastic strain
      Ceps      = data(13); // equivalent plastic strain
      Cdamage   = data(14); // damage variable controlling compressive softening     

      // Copy converged state values into trial values
      Tstrain   = Cstrain;   // strain
      Tnlstrain = Cnlstrain; // non-local strain
      Tstress   = Cstress;   // stress
      Ttangent  = Ctangent;  // material stiffness (tangent)
      Tpstrain  = Cpstrain;  // total plastic strain
      Teps      = Ceps;      // equivalent plastic strain
      Tdamage   = Cdamage;   // damage variable controlling compressive softening
   }
    
   return res;
}

void Steel01::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "Steel01 tag: " << this->getTag() << endln;
    s << "  fy: " << fy << " ";
    s << "  E0: " << E0 << " ";
    s << "   b: " << b << " ";
    s << "  a1: " << a1 << " ";
    s << "  a2: " << a2 << " ";
    s << "  a3: " << a3 << " ";
    s << "  a4: " << a4 << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"Steel01\", ";
	s << "\"E\": " << E0 << ", ";
	s << "\"fy\": " << fy << ", ";
    s << "\"b\": " << b << ", ";
    s << "\"a1\": " << a1 << ", ";
    s << "\"a2\": " << a2 << ", ";
    s << "\"a3\": " << a3 << ", ";
    s << "\"a4\": " << a4 << "}";
  }
  
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
Steel01::setParameter(const char **argv, int argc, Parameter &param)
{
        return 0;
}



int
Steel01::updateParameter(int parameterID, Information &info)
{
	return 0;
}




int
Steel01::activateParameter(int passedParameterID)
{
	return 0;
}



double
Steel01::getStressSensitivity(int gradIndex, bool conditional)
{
	return 0.0;
}




double
Steel01::getInitialTangentSensitivity(int gradIndex)
{
        return 0.0;
}


int
Steel01::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////

