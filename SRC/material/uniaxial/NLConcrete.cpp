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

// $Revision: 1.1 $
// $Date: 2009-05-29 19:10:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NLConcrete.cpp,v $

// Written: Maha Kenawy
// Created: August 2020

// refer to paper:
//


#include <NLConcrete.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>
#include <math.h>
#include <float.h>
#include <MatrixOperations.h>



NLConcrete::NLConcrete(int mTag, double mE, double mfc,
				     double mec0, double mEd)
:UniaxialMaterial(mTag,MAT_TAG_NLConcrete),
 E(mE), fc(mfc), ec0(mec0), Ed(mEd)
{
	// Initialize variables
    this->revertToStart();
}

NLConcrete::~NLConcrete()
{
}

int NLConcrete::setNLStrain(double nlstrain)
{
    Tnlstrain = nlstrain;
    return 0;
}

int 
NLConcrete::setTrialStrain (double strain, double strainRate)
{
	// material properties
	TStrain = strain;
	double dStrain = TStrain - CStrain;
	double dnlstrain = Tnlstrain - Cnlstrain;
	double fy = 0.5*fc;
	double Tystress = fy;
	double k0 = fabs(ec0 - fc/E);
	double H = (fabs(fc) - fabs(fy))/k0;
	double Mod = H;
	// linear damage constants
	double ef = ec0 - fc/Ed;
	double kf = fabs(ef - fc/E);
	double alphalin = 1.0;
	double H2 = 0.15*E;

    // compute elastic trial stress
	double deltaLambda = 0.0;
	TPlasticStrain = CPlasticStrain;
	TBackStress = CBackStress;
	TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
	TD = CD;
	Tnlk = Cnlk;
	Tkd2 = Ckd2;
	double dDoverdk = 0.0;
	double m = 1.0;
	double efs = ef - k0;
	int elastic = 1;

    // beginning of plasticity
	double TrialStresstmp = E*(TStrain - CPlasticStrain);
	double TStresstmp = TrialStresstmp;

	// no tension is accommodated
	 if (TrialStresstmp > 0.0)
	{
		TStress = 0.0;
		TTangent = 1.0e-10;
	//        TTangent = E;
	}

	else {
		// trial plastic stress
		 if (TAccumulatedPlasticStrain <= k0)
		 {
			Mod = H;
			Tystress = fy + Mod*TAccumulatedPlasticStrain;
		 }
		 else {
			Mod = H2;
			Tystress = fc + Mod*(TAccumulatedPlasticStrain - k0);
		 }

    // Compute yield criterion
    double f = fabs(TrialStresstmp) - Tystress;

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
		// Set trial tangent
		TTangent = (1 - TD)*E;
        elastic = 1;
    }

    // Plastic step ... perform return mapping algorithm ---
    else {
        elastic = 0;
        deltaLambda = f/(E+Mod);

        // Find sign of stress
        int sign = (TStresstmp < 0) ? -1 : 1;

        TPlasticStrain = CPlasticStrain +deltaLambda*sign;
        TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;
        TStresstmp = TrialStresstmp - E*deltaLambda*sign;
        double estraininc = dStrain - deltaLambda*sign;

	  // calculate the damage
      // calculate the nonlocal damage variables kd1 and kd2
	  double dnlk = std::max(fabs(dnlstrain) - fabs(estraininc),0.0);
	  Tnlk = Cnlk + dnlk;
	  //double dnlk = std::max(deltaLambda + m*fabs(dnlstrain) - m*fabs(dStrain),0.0);
	  //Tnlk = TAccumulatedPlasticStrain + m*fabs(Tnlstrain) - m*fabs(TStrain);
	  if (Tnlk > k0)
	  {
		  double kd1 = Tnlk - k0;
		  TD = (E*kd1*(efs*H2 + fc))/((fc + H2*kd1)*(E*efs - fc));
		  dDoverdk = fc*E*(efs*H2 + fc)/(E*efs - fc)/pow((fc + H2*kd1),2);
	  }

	  if (TD > 1.0)
		  TD = 1.0;
	  if (fabs(TD - CD) < DBL_EPSILON)
		  dDoverdk = 0.0;

      TTangent = ((1 - TD)*Mod/(E+Mod) - TStresstmp*dDoverdk*-1/(E+Mod))*E;
      if (TD > 0.999)
    	  TTangent = 1.0e-10;
    }
    TStress = (1 - TD)*TStresstmp;
        }
    return 0;
}

double 
NLConcrete::getStress(void)
{
    return TStress;
}

double 
NLConcrete::getTangent(void)
{
    return TTangent;
}

double 
NLConcrete::getStrain(void)
{
    return TStrain;
}

int 
NLConcrete::commitState(void)
{
    // Commit trial history variables
    CPlasticStrain = TPlasticStrain;
    CBackStress = TBackStress;
    CAccumulatedPlasticStrain = TAccumulatedPlasticStrain;

    CStrain = TStrain;		// Committed strain
    Cnlstrain = Tnlstrain;
    CStress = TStress;		// Committed stress
    CTangent = TTangent;
    CD = TD;
    Cnlk = Tnlk;
    Ckd2 = Tkd2;
    
    return 0;
}

int 
NLConcrete::revertToLastCommit(void)
{
  return 0;
}

int 
NLConcrete::revertToStart(void)
{
    // Reset committed history variables
    CPlasticStrain = 0.0;
    CBackStress = 0.0;
    CAccumulatedPlasticStrain = 0.0;

    // Reset committed history variables
    TPlasticStrain = 0.0;
    TBackStress = 0.0;
    TAccumulatedPlasticStrain = 0.0;

	// Initialize state variables
	TStrain = 0.0;
    Tnlstrain = 0.0;
	TStress = 0.0;
	TTangent = E;
	TD = 0.0;
	Tnlk = 0.0;
	Tkd2 = 0.0;

	// Initialize committed variables
	CStrain = 0.0;
	Cnlstrain = 0.0;
	CStress = 0.0;
	CTangent = E;
	CD = 0.0;
	Cnlk = 0.0;
	Ckd2 = 0.0;
	maxstrain = 0.0;

    return 0;
}

UniaxialMaterial *
NLConcrete::getCopy(void)
{
    NLConcrete *theCopy =
	new NLConcrete(this->getTag(), E, fc, ec0, Ed);

    // Copy committed history variables
    theCopy->CPlasticStrain = CPlasticStrain;
    theCopy->CBackStress = CBackStress;
    theCopy->CAccumulatedPlasticStrain = CAccumulatedPlasticStrain;

    // Copy trial history variables
    theCopy->TPlasticStrain = TPlasticStrain;
    theCopy->TBackStress = TBackStress;
    theCopy->TAccumulatedPlasticStrain = TAccumulatedPlasticStrain;

    // Copy trial state variables
    theCopy->TStrain = TStrain;
    theCopy->Tnlstrain = Tnlstrain;
    theCopy->TStress = TStress;
    theCopy->TTangent = TTangent;
    theCopy->CStrain = CStrain;
    theCopy->Cnlstrain = Cnlstrain;
    theCopy->CStress = CStress;
    theCopy->CTangent = CTangent;
    theCopy->CD = CD;
    theCopy->TD = TD;
    theCopy->Cnlk = Cnlk;
    theCopy->Tnlk = Tnlk;
    theCopy->Ckd2 = Ckd2;
    theCopy->Tkd2 = Tkd2;
    theCopy->maxstrain = maxstrain;

    return theCopy;
}

int 
NLConcrete::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(23);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = fc;
  data(3) = ec0;
  data(4) = Ed;
  data(5) = CPlasticStrain;
  data(6) = CBackStress;
  data(7) = CAccumulatedPlasticStrain;
  data(8) = TStrain;
  data(9) = Tnlstrain;
  data(10) = TStress;
  data(11) = TTangent;
  data(12) = CStrain;
  data(13) = Cnlstrain;
  data(14) = CStress;
  data(15) = CTangent;
  data(16) = TD;
  data(17) = CD;
  data(18) = Tnlk;
  data(19) = Cnlk;
  data(20) = Tkd2;
  data(21) = Ckd2;
  data(22) = maxstrain;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "NLConcrete::sendSelf() - failed to send data\n";

  return res;
}

int 
NLConcrete::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(23);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "NLConcrete::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    fc = data(2);
    ec0 = data(3);
    Ed = data(4);
    CPlasticStrain = data(5);
    CBackStress = data(6);
    CAccumulatedPlasticStrain = data(7);
    TStrain = data(8);
    Tnlstrain = data(9);
    TStress = data(10);
    TTangent = data(11);
    CStrain = data(12);
    Cnlstrain = data(13);
    CStress = data(14);
    CTangent = data(15);
    TD = data(16);
    CD = data(17);
    Tnlk = data(18);
    Cnlk = data(19);
    Tkd2 = data(20);
    Ckd2 = data(21);
    maxstrain = data(22);
  }
    
  return res;
}

void 
NLConcrete::Print(OPS_Stream &s, int flag)
{
    s << "NLConcrete, tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
    s << "  fc: " << fc << endln;
    s << "  ec0: " << ec0 << endln;
    s << "  Ed: " << Ed << endln;
   
}
