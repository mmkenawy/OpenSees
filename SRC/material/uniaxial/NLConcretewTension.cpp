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

// $Revision: 2.0 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NLConcretewTension.cpp,v $

// Author: Maha Kenawy
// Created: August 2020

// refer to paper:
//


#include <NLConcretewTension.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>
#include <math.h>
#include <float.h>
#include <MatrixOperations.h>


NLConcretewTension::NLConcretewTension(int pTag, double pE, double pfc,
				     double pec0, double pEd, double pft, double peft)
:UniaxialMaterial(pTag,MAT_TAG_NLConcretewTension),
 E(pE), fc(pfc), ec0(pec0), Ed(pEd), ft(pft), eft(peft)
{
	// Initialize variables
    this->revertToStart();
}

NLConcretewTension::~NLConcretewTension()
{
}

int NLConcretewTension::setNLStrain(double nlstrain)
{
    Tnlstrain = nlstrain;

    return 0;
}

int 
NLConcretewTension::setTrialStrain (double strain, double strainRate)
{
	// Define some material parameters
	// strain increments
	TStrain = strain;
	double dStrain = TStrain - CStrain;
	double dnlstrain = Tnlstrain - Cnlstrain;
	// compressive properties
	double fy = 0.5*fc;
	double sigyc = fy;
	double k0 = fabs(ec0 - fc/E);
	double H = (fabs(fc) - fabs(fy))/k0;
	double Hc = H;
	double H2 = 0.15*E;
	double ef = ec0 - fc/Ed;
	double efs = ef - k0;

	// tensile properties
	//double ft = 0.1*fc;
	double sigyt = ft;
	//double eft = 0.003;
	double Ht = 0.15*E;

	// crack closing modulus
	double H3 = 0.05*E;

    // ------ Elastic trial + crack-closing -------
	double deltaLambda = 0.0;
	int elastic = 1;
	TPlasticStrain = CPlasticStrain;
    Tefstress = Cefstress;
    Tkk = Ckk;
	Tkc = Ckc;
	Tkt = Ckt;
	TDc = CDc;
	TDt = CDt;
	Tkdc1 = Ckdc1;
	Tkdt1 = Ckdt1;
	TDam = CDam;
	double dDdkdc = 0.0;
	double dDdkdt = 0.0;

	double TrialStresstmp = E*(TStrain - TPlasticStrain);

	double yalpha = 0.0; // material under compression
	if (TrialStresstmp > 0.0) // material under tension
	{
			yalpha = 1.0;
     }

//  crack closing model
	double crackinc = 0.0;
	double dcrack = E/(E + H3);
	if (TrialStresstmp < 0.0 && dStrain < 0.0 && Tkt + Tkk > 0.0) { // if the material is in crack-closing phase
			crackinc = dcrack*dStrain;
			TPlasticStrain += crackinc;
			Tkk = Ckk + crackinc;
			TrialStresstmp -= E*crackinc;
				  }

// compute the value of the current yield stresses in tension and compression
	// current tensile yield stress
     sigyt = ft + Ht*Ckt;

     // current compressive yield stress
	 if (Ckc <= k0)
	 {
		sigyc = fy + H*Ckc;
		Hc = H;
	 }
	 else {
		sigyc = fc + H2*(Ckc - k0);
		Hc = H2;
	 }

	 // current effective material modulus
	 double effmod = pow((1 - yalpha),2)*Hc + pow(yalpha,2)*Ht;

    // Compute yield criterion
    double f = fabs(TrialStresstmp) - (1-yalpha)*sigyc - yalpha*sigyt;
    int sign = (TrialStresstmp < 0) ? -1 : 1;

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
        Tefstress = TrialStresstmp;
    }

    // ------- Plastic step ... perform return mapping algorithm ---
    else {
    	elastic = 0;
        deltaLambda = f/(E+effmod);

	    TPlasticStrain += deltaLambda*sign;
	    Tefstress = TrialStresstmp - E*deltaLambda*sign;

	    // internal damage variables
        Tkc = Ckc + (1-yalpha)*deltaLambda;
        Tkt = Ckt + yalpha*deltaLambda;
        //elastic strain increment
        double estraininc = dStrain - deltaLambda*sign;

        // evaluate damage in tension
		double dnlkt = yalpha*std::max(fabs(dnlstrain) - fabs(estraininc),0.0);
		Tkdt1 = Ckdt1 + dnlkt;
		if (Tkdt1 > 0.0)
        {
		  TDt = std::min((E*Tkdt1*(eft*Ht + ft))/((ft + Ht*Tkdt1)*(E*eft - ft)),1.0);
		  dDdkdt = ft*E*(eft*Ht + ft)/(E*eft - ft)/pow((ft + Ht*Tkdt1),2);

          }
		// evaluate damage in compression
		double dnlkc = (1-yalpha)*std::max(fabs(dnlstrain) - fabs(estraininc),0.0);
		Tkdc1 = Ckdc1 + dnlkc;
		if (Tkdc1 > k0)
		{
		  TDc = std::min((E*(Tkdc1 - k0)*(efs*H2 + fc))/((fc + H2*(Tkdc1 - k0))*(E*efs - fc)),1.0);
		  dDdkdc = fc*E*(efs*H2 + fc)/(E*efs - fc)/pow((fc + H2*(Tkdc1 - k0)),2);
		  }
    	}

    // total material damage at this step
    TDam = (1-yalpha)*TDc + yalpha*TDt;

    // derivative of damage (needed for tangent computation)
    double dDeff = pow((1-yalpha),2)*dDdkdc + pow(yalpha,2)*dDdkdt;

    if (TDam > 0.99999 || TDam - CDam < DBL_EPSILON)
       dDeff = 0.0;

    // nominal stress
    TStress = (1 - TDam)*Tefstress;

    // compute material tangent
    if (elastic == 1) // if material is in elastic range
    {
        if (Tefstress < 0.0 && dStrain < 0.0 && Tkt+ Tkk > 0.0) // if material is in crack-closing phase
				TTangent = (1-TDam)*E*H3/(E + H3);
		else
				TTangent = (1-TDam)*E;
    }
    else // if material is in plasticity/damage range
    {
        TTangent = ((1 - TDam)*effmod - Tefstress*dDeff*sign)*E/(E+effmod);
        //TTangent = (1 - TDam)*E;
    }

    return 0;
}

double 
NLConcretewTension::getStress(void)
{
    return TStress;
}

double 
NLConcretewTension::getTangent(void)
{
    return TTangent;
}

double 
NLConcretewTension::getStrain(void)
{
    return TStrain;
}

int 
NLConcretewTension::commitState(void)
{
	// Commit trial history variables
	CPlasticStrain = TPlasticStrain;
	CDam = TDam;
	Ckk = Tkk;
	Ckc = Tkc;
	Ckt = Tkt;
	CStrain = TStrain;		// Committed strain
	Cnlstrain = Tnlstrain;
	CStress = TStress;		// Committed stress
	CTangent = TTangent;
	CDc = TDc;
	CDt = TDt;
	Cefstress = Tefstress;
	Ckdc1 = Tkdc1;
	Ckdt1 = Tkdt1;
    
	return 0;
}

int 
NLConcretewTension::revertToLastCommit(void)
{
  return 0;
}

int 
NLConcretewTension::revertToStart(void)
{
    // Reset committed history variables
    CPlasticStrain = 0.0;
    CDam = 0.0;
    Ckk = 0.0;
    Ckc = 0.0;
    Ckt = 0.0;

    // Reset committed history variables
    TPlasticStrain = 0.0;
    TDam = 0.0;
    Tkk = 0.0;
    Tkc = 0.0;
    Tkt = 0.0;

	// Initialize state variables
	TStrain = 0.0;
    Tnlstrain = 0.0;
	TStress = 0.0;
    Tefstress = 0.0;
	TTangent = E;
    TDc = 0.0;
    TDt = 0.0;

	// Initialize committed variables
	CStrain = 0.0;
    Cnlstrain = 0.0;
	CStress = 0.0;
    Cefstress = 0.0;
	CTangent = E;
	CDc = 0.0;
	CDt = 0.0;
	Ckdc1 = 0.0;
	Ckdt1 = 0.0;

    return 0;
}

UniaxialMaterial *
NLConcretewTension::getCopy(void)
{
    NLConcretewTension *theCopy =
	new NLConcretewTension(this->getTag(), E, fc, ec0, Ed, ft, eft);

    // Copy committed history variables
    theCopy->CPlasticStrain = CPlasticStrain;
    theCopy->Ckc = Ckc;
    theCopy->Ckt = Ckt;

    // Copy trial history variables
    theCopy->TPlasticStrain = TPlasticStrain;
    theCopy->Tkc = Tkc;
    theCopy->Tkt = Tkt;

    // Copy trial state variables
    theCopy->TStrain = TStrain;
    theCopy->Tnlstrain = Tnlstrain;
    theCopy->TStress = TStress;
    theCopy->TTangent = TTangent;
    theCopy->CStrain = CStrain;
    theCopy->Cnlstrain = Cnlstrain;
    theCopy->CStress = CStress;
    theCopy->CTangent = CTangent;
    theCopy->CDc = CDc;
    theCopy->TDc = TDc;
    theCopy->CDt = CDt;
    theCopy->TDt = TDt;
    theCopy->Tefstress = Tefstress;
    theCopy->Cefstress = Cefstress;
    theCopy->Tkdc1 = Tkdc1;
    theCopy->Tkdt1 = Tkdt1;
    theCopy->TDam = TDam;
    theCopy->CDam = CDam;
    theCopy->Tkk = Tkk;
    theCopy->Ckk = Ckk;
    return theCopy;
}

int 
NLConcretewTension::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(31);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = fc;
  data(3) = ec0;
  data(4) = Ed;
  data(5) = ft;
  data(6) = eft;
  data(7) = CPlasticStrain;
  data(8) = Ckc;
  data(9) = Ckt;
  data(10) = TStrain;
  data(11) = Tnlstrain;
  data(12) = TStress;
  data(13) = TTangent;
  data(14) = CStrain;
  data(15) = Cnlstrain;
  data(16) = CStress;
  data(17) = CTangent;
  data(18) = CDc;
  data(19) = CDt;
  data(20) = Tkc;
  data(21) = Tkt;
  data(22) = TDc;
  data(23) = TDt;
  data(24) = Cefstress;
  data(25) = Tefstress;
  data(26) = Ckdc1;
  data(27) = Ckdt1;
  data(28) = CDam;
  data(29) = TDam;
  data(30) = Ckk;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "NLConcretewTension::sendSelf() - failed to send data\n";

  return res;
}

int 
NLConcretewTension::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(31);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "NLConcretewTension::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    fc = data(2);
    ec0 = data(3);
    Ed = data(4);
    ft = data(5);
    eft = data(6);
    CPlasticStrain = data(7);
    Ckc = data(8);
    Ckt = data(9);
    TStrain = data(10);
    Tnlstrain = data(11);
    TStress = data(12);
    TTangent = data(13);
    CStrain = data(14);
    Cnlstrain = data(15);
    CStress = data(16);
    CTangent = data(17);
    CDc = data(18);
    CDt = data(19);
    Tkc = data(20);
    Tkt = data(21);
    TDc = data(22);
    TDt = data(23);
    Cefstress = data(24);
    Tefstress = data(25);
    Ckdc1 = data(26);
    Ckdt1 = data(27);
    CDam = data(28);
    TDam = data(29);
    Ckk = data(30);
  }
    
  return res;
}

void 
NLConcretewTension::Print(OPS_Stream &s, int flag)
{
    s << "NLConcretewTension, tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
    s << "  fc: " << fc << endln;
    s << "  ec0: " << ec0 << endln;
    s << "  Ed: " << Ed << endln;
    s << "  ft: " << ft << endln;
    s << "  eft: " << eft << endln;
   
}
