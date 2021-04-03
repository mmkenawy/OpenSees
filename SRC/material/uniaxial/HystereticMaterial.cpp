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
                                                                        
// $Revision: 1.19 $
// $Date: 2008-12-18 23:40:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticMaterial.cpp,v $

// Written: MHS
// Created: July 2000
//
// Description: This file contains the implementation of 
// HystereticMaterial.  HystereticMaterial is
// a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.
#include <stdlib.h>
#include <HystereticMaterial.h>
#include <OPS_Globals.h>
#include <math.h>
#include <float.h>
#include <Channel.h>

#include <elementAPI.h>

void *
OPS_HystereticMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 18 && numArgs != 17 && numArgs != 14 && numArgs != 13) {
    opserr << "Want: uniaxialMaterial Hysteretic tag? mom1p? rot1p? mom2p? rot2p? <mom3p? rot3p?> "
	   << "\nmom1n? rot1n? mom2n? rot2n? <mom3n? rot3n?> pinchX? pinchY? damfc1? damfc2? <beta?>";
    return 0;
  }
  
  int iData[1];
  double dData[17];
  for (int i=0; i<17; i++) 
    dData[i] = 0.0;

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Hysteretic" << endln;
    return 0;
  }

  numData = numArgs-1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial Hysteretic " << iData[0] << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  if (numData > 13) 
    theMaterial = new HystereticMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
					 dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
  					 dData[13], dData[14], dData[15], dData[16]);
  else
    theMaterial = new HystereticMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
					 dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Hysteretic\n";
    return 0;
  }

  return theMaterial;
}



HystereticMaterial::HystereticMaterial(int tag,
				       double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
				       double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
				       double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom2p(m2p), rot2p(r2p), mom3p(m3p), rot3p(r3p),
mom1n(m1n), rot1n(r1n), mom2n(m2n), rot2n(r2n), mom3n(m3n), rot3n(r3n)
{
  bool error = false;
  // Positive backbone parameters
  if (rot1p <= 0.0)
    error = true;
  
  if (rot2p <= rot1p)
    error = true;
  
  if (rot3p <= rot2p)
    error = true;
  
  // Negative backbone parameters
  if (rot1n >= 0.0)
    error = true;
  
  if (rot2n >= rot1n)
    error = true;
  
  if (rot3n >= rot2n)
    error = true;
  
  if (error) {
    opserr << "HystereticMaterial::HystereticMaterial -- input backbone is not unique (one-to-one)\n";
    exit(-1);
  }		
  
  energyA = 0.5 * (rot1p*mom1p + (rot2p-rot1p)*(mom2p+mom1p) + (rot3p-rot2p)*(mom3p+mom2p) +
		   rot1n*mom1n + (rot2n-rot1n)*(mom2n+mom1n) + (rot3n-rot2n)*(mom3n+mom2n));
  
  // Set envelope slopes
  this->setEnvelope();
  
  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
  
}

HystereticMaterial::HystereticMaterial(int tag,
			double m1p, double r1p, double m2p, double r2p,
			double m1n, double r1n, double m2n, double r2n,
			double px, double py, double d1, double d2, double b):
UniaxialMaterial(tag, MAT_TAG_Hysteretic),
pinchX(px), pinchY(py), damfc1(d1), damfc2(d2), beta(b),
mom1p(m1p), rot1p(r1p), mom3p(m2p), rot3p(r2p),
mom1n(m1n), rot1n(r1n), mom3n(m2n), rot3n(r2n)
{
	bool error = false;
	
	// Positive backbone parameters
	if (rot1p <= 0.0)
		error = true;

	if (rot3p <= rot1p)
		error = true;

	// Negative backbone parameters
	if (rot1n >= 0.0)
		error = true;

	if (rot3n >= rot1n)
		error = true;

	if (error) {
	  opserr << "HystereticMaterial::HystereticMaterial -- input backbone is not unique (one-to-one)\n";
	  exit(-1);
	}

				      

	energyA = 0.5 * (rot1p*mom1p + (rot3p-rot1p)*(mom3p+mom1p) +
		rot1n*mom1n + (rot3n-rot1n)*(mom3n+mom1n));

	mom2p = 0.5*(mom1p+mom3p);
	mom2n = 0.5*(mom1n+mom3n);

	rot2p = 0.5*(rot1p+rot3p);
	rot2n = 0.5*(rot1n+rot3n);

	// Set envelope slopes
	this->setEnvelope();

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

HystereticMaterial::HystereticMaterial():
UniaxialMaterial(0, MAT_TAG_Hysteretic),
pinchX(0.0), pinchY(0.0), damfc1(0.0), damfc2(0.0), beta(0.0),
mom1p(0.0), rot1p(0.0), mom2p(0.0), rot2p(0.0), mom3p(0.0), rot3p(0.0),
mom1n(0.0), rot1n(0.0), mom2n(0.0), rot2n(0.0), mom3n(0.0), rot3n(0.0)
{

}

HystereticMaterial::~HystereticMaterial()
{
	// Nothing to do
}

int
HystereticMaterial::setTrialStrain(double strain, double strainRate)
{
  if (TloadIndicator == 0 && strain == 0.0)
    return 0;

  // Re-set trial values to their original converged values from previous step
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TenergyD = CenergyD;
  TrotPu = CrotPu;
  TrotNu = CrotNu;
  TnlstrainMax = CnlstrainMax;
  TnlstrainMin = CnlstrainMin;

  // Set the trial strain value in the material, and compute the strain increment
  Tstrain = strain;
  double dStrain = Tstrain - Cstrain;

  // Do not update the un-damanged stress state if the strain increment is close to zero
  if (fabs(dStrain) < DBL_EPSILON) {
    Tustress = Custress;
  } else {
  
    // Re-set load indicator
    TloadIndicator = CloadIndicator;

    // If this is the first step, set the load indicator to:
    // 1: if strain increment is positive or zero
    // 2: if strain increment is negative
    if (TloadIndicator == 0) TloadIndicator = (dStrain < 0.0) ? 2 : 1;

    // Update the undamaged stress:
    // Check whether the trial strain exceeds certain min/max strain bounds
    if (Tstrain >= CrotMax) {
      // if strain exceeds max value, assume we are on the (positive) +envelope
      TrotMax = Tstrain;
      // compute tangent and stress based upon a given strain position on the +envelope
      Tutangent = posEnvlpTangent(Tstrain);
      Tustress = posEnvlpStress(Tstrain);
      // set the load indicator to 1: positive loading on the +envelope
      TloadIndicator=1;
    }
    else if (Tstrain <= CrotMin) {
      // if strain is below min value, assume we are on the (negative) -envelope
      TrotMin = Tstrain;
      // compute tangent and stress based upon a given strain position on the -envelope
      Tutangent = negEnvlpTangent(Tstrain);
      Tustress = negEnvlpStress(Tstrain);
      // set the load indicator to 2: negative loading on the -envelope
      TloadIndicator=2;
    }
    else {
      // otherwise, assume we are somewhere inside of the quasi-"elastic" range
      if (dStrain < 0.0)
	// negative elastic (un)loading
	negativeIncrement(dStrain);
      else if (dStrain > 0.0)
	// positive elastic loading
	positiveIncrement(dStrain);
    }
  }

  // independently set the min/max trial non-local strain:
  if (Tnlstrain >= CnlstrainMax) TnlstrainMax = Tnlstrain;
  if (Tnlstrain <= CnlstrainMin) TnlstrainMin = Tnlstrain;

  // Re-scale the stress by the damage factor in tension/compression:
  if (Tustress < 0.0) {
    // compute current damage in compression
    double m = 1.5;
    double rot0n = -mom2n/E3n + rot2n;
    double nlrot = m*TnlstrainMin + (1.0-m)*TrotMin;
    double dam = 1.0 - (fabs(rot0n)-fabs(nlrot))/(fabs(rot0n)-fabs(rot2n));
    if (dam > 0.0)
          double dum = 0.0;
    if (dam < 0.0)
      dam = 0.0;
    if (dam > 0.8)
      dam = 0.8;
    // apply negative damage scaling
    Tstress = (1.0-dam)*Tustress;
    Ttangent = (1.0-dam)*Tutangent;
  } else {
    // compute current damage in tension
    double m = 1.5;
    double rot0p = -mom2p/E3p + rot2p;
    double nlrot = m*TnlstrainMax + (1.0-m)*TrotMax;
    double dam = 1.0 - (fabs(rot0p)-fabs(nlrot))/(fabs(rot0p)-fabs(rot2p));
    if (dam > 0.0)
      double dum = 0.0;
    if (dam < 0.0)
      dam = 0.0;
    if (dam > 0.8)
      dam = 0.8;
    // apply positive damage scaling
    Tstress = (1.0-dam)*Tustress;
    Ttangent = (1.0-dam)*Tutangent;
  }

  // time-integrate: (trial internal energy) = time integral of { (stress) * (strain rate) }
  TenergyD = CenergyD + 0.5*(Custress+Tustress)*dStrain;

  //  if (this->getTag() == 40)
  //    opserr << "setTrial: " << Tstrain << " " << Tutangent << " " << Tstress << endln;
  
  return 0;
}

int
HystereticMaterial::setNLStrain(double nlstrain)
{
    Tnlstrain = nlstrain;
    return 0;
}

double
HystereticMaterial::getStrain(void)
{
	return Tstrain;
}

double
HystereticMaterial::getStress(void)
{
	return Tstress;
}

double
HystereticMaterial::getTangent(void)
{
  return Ttangent;
}

void
HystereticMaterial::positiveIncrement(double dStrain)
{
	double kn = pow(CrotMin/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;
	
	// if the loading indicator shows a change from negative loading to positive loading
	// (as indicated by the positive sign of the strain increment)
	// then flip the load indicator to indicate positive loading
	if (TloadIndicator == 2) {
		TloadIndicator = 1;
		// if the stress is currently less than zero, compute new TrotNu (plastic strain)
		if (Custress <= 0.0) {
		  // solve for TrotNu => Custress = Eun*kn*(Cstrain - TrotNu)
		  // TrotNu is essentially the current reference 'plastic strain' value, e.g.
		  // s = E*(e - e^p)
			TrotNu = Cstrain - Custress/(Eun*kn);
			double energy = CenergyD - 0.5*Custress/(Eun*kn)*Custress;
			double damfc = 0.0;
			if (CrotMin < rot1n) {
				damfc = damfc2*energy/energyA;
				damfc += damfc1*(CrotMin-rot1n)/rot1n;
			}
			// re-load the max strain (without damage, this should be unaltered)

			TrotMax = CrotMax*(1.0+damfc);
		}
	}

	// set the load indicator to match the positive direction of the strain increment
	TloadIndicator = 1;

	// threshhold the maximum strain (possibly for improving convergence or managing overflow?)
  if (TrotMax > POS_INF_STRAIN)
    TrotMax = POS_INF_STRAIN;

  // threshhold the trial max strain to the lower bound of rot1p
	TrotMax = (TrotMax > rot1p) ? TrotMax : rot1p;

	// evaluate the positive envelope stress at the trial max strain value
	double maxmom = posEnvlpStress(TrotMax);
	// evaluate the negative strain limit: the min allowable value for the min strain value
	double rotlim = negEnvlpRotlim(CrotMin);
	// evaluate the relative strain as the larger of the negative strain limit or TrotNu
	// (if the min allowable strain is numerical -infinity, rotrel will always be TrotNu)
	double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;

	// rotrel = TrotNu;
	// if (negEnvlpStress(CrotMin) >= 0.0)
	//    rotrel = rotlim;
	
	//	double rotmp1 = rotrel + pinchY*(TrotMax-rotrel);

	// Compute the pinching strain value (during unloading?) based upon the projected envelope stress
	// solve for rotmp2 => (Eup*kp)*(TrotMax - rotmp2) = (1.0-pinchY)*maxmom
	double rotmp2 = TrotMax - (1.0-pinchY)*maxmom/(Eup*kp);
	//double rotmp2 = TrotMax-(1-pinchY)*maxmom/Eup;
	//	double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
	// if pinchX in the range of 0 to 1 and rotrel=TrotNu, rotch is the pinching strain during loading?
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;                   // changed on 7/11/2006

	double tmpmo1;
	double tmpmo2;

	// if the trial strain is less than the current plastic strain (negative stress regime):
	if (Tstrain < TrotNu) {
		Tutangent = Eun*kn;
		Tustress = Custress + Tutangent*dStrain;
		if (Tustress >= 0.0) {
			Tustress = 0.0;
			Tutangent = Eun*1.0e-9;
		}
	}
	// if the trial strain exceeds the current plastic strain (positive stress regime)
	// and we are on the first branch of the pinching curve
	else if (Tstrain >= TrotNu && Tstrain < rotch) {
		if (Tstrain <= rotrel) {
			Tustress = 0.0;
			Tutangent = Eup*1.0e-9;
		}
		// increment the stress along the first branch of the pinching curve
		else {
			Tutangent = maxmom*pinchY/(rotch-rotrel);
			tmpmo1 = Custress + Eup*kp*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Tutangent;
			if (tmpmo1 < tmpmo2) {
				Tustress = tmpmo1;
				Tutangent = Eup*kp;
			}
			else
				Tustress = tmpmo2;
		}
	}
	// if we are on the second branch of the pinching curve, increment stress accordingly
	else {
		Tutangent = (1.0-pinchY)*maxmom/(TrotMax-rotch);
		tmpmo1 = Custress + Eup*kp*dStrain;
		tmpmo2 = pinchY*maxmom + (Tstrain-rotch)*Tutangent;
		if (tmpmo1 < tmpmo2) {
			Tustress = tmpmo1;
			Tutangent = Eup*kp;
		}
		else
			Tustress = tmpmo2;
	}
}

void
HystereticMaterial::negativeIncrement(double dStrain)
{
	double kn = pow(CrotMin/rot1n,beta);
	kn = (kn < 1.0) ? 1.0 : 1.0/kn;
	double kp = pow(CrotMax/rot1p,beta);
	kp = (kp < 1.0) ? 1.0 : 1.0/kp;

	if (TloadIndicator == 1) {
		TloadIndicator = 2;
		if (Custress >= 0.0) {
			TrotPu = Cstrain - Custress/(Eup*kp);
			double energy = CenergyD - 0.5*Custress/(Eup*kp)*Custress;
			double damfc = 0.0;
			if (CrotMax > rot1p) {
				damfc = damfc2*energy/energyA;
				damfc += damfc1*(CrotMax-rot1p)/rot1p;
			}

			TrotMin = CrotMin*(1.0+damfc);
		}
	}

  TloadIndicator = 2;

  if (TrotMin < NEG_INF_STRAIN)
    TrotMin = NEG_INF_STRAIN;

	TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;

	double minmom = negEnvlpStress(TrotMin);
	double rotlim = posEnvlpRotlim(CrotMax);
	double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;

	//rotrel = TrotPu;
	//if (posEnvlpStress(CrotMax) <= 0.0)
	//  rotrel = rotlim;

	//double rotmp1 = rotrel + pinchY*(TrotMin-rotrel);
	double rotmp2 = TrotMin - (1.0-pinchY)*minmom/(Eun*kn);
	//double rotmp2 = TrotMin-(1-pinchY)*minmom/Eun;	
	//double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
	double rotch = rotrel + (rotmp2-rotrel)*pinchX;                   // changed on 7/11/2006

	double tmpmo1;
	double tmpmo2;

	if (Tstrain > TrotPu) {
		Tutangent = Eup*kp;
		Tustress = Custress + Tutangent*dStrain;
		if (Tustress <= 0.0) {
			Tustress = 0.0;
			Tutangent = Eup*1.0e-9;
		}
	}

	else if (Tstrain <= TrotPu && Tstrain > rotch) {
		if (Tstrain >= rotrel) {
			Tustress = 0.0;
			Tutangent = Eun*1.0e-9;
		}
		else {
			Tutangent = minmom*pinchY/(rotch-rotrel);
			tmpmo1 = Custress + Eun*kn*dStrain;
			tmpmo2 = (Tstrain-rotrel)*Tutangent;
			if (tmpmo1 > tmpmo2) {
				Tustress = tmpmo1;
				Tutangent = Eun*kn;
			}
			else
				Tustress = tmpmo2;
		}
	}

	else {
		Tutangent = (1.0-pinchY)*minmom/(TrotMin-rotch);
		tmpmo1 = Custress + Eun*kn*dStrain;
		tmpmo2 = pinchY*minmom + (Tstrain-rotch)*Tutangent;
		if (tmpmo1 > tmpmo2) {
			Tustress = tmpmo1;
			Tutangent = Eun*kn;
		}
		else
			Tustress = tmpmo2;
	}
}

int
HystereticMaterial::commitState(void)
{
	CrotMax = TrotMax;
	CrotMin = TrotMin;
	CrotPu = TrotPu;
	CrotNu = TrotNu;
	CenergyD = TenergyD;
	CloadIndicator = TloadIndicator;

	Cstress = Tstress;
	Custress = Tustress;
	Cstrain = Tstrain;
	Cnlstrain = Tnlstrain;
	CnlstrainMax = TnlstrainMax;
	CnlstrainMin = TnlstrainMin;
	return 0;
}

int
HystereticMaterial::revertToLastCommit(void)
{
	TrotMax = CrotMax;
	TrotMin = CrotMin;
	TrotPu = CrotPu;
	TrotNu = CrotNu;
	TenergyD = CenergyD;
	TloadIndicator = CloadIndicator;

	Tstress = Cstress;
	Tustress = Custress;
	Tstrain = Cstrain;
	Tnlstrain = Cnlstrain;
	TnlstrainMax = CnlstrainMax;
	TnlstrainMin = CnlstrainMin;

	return 0;
}

int
HystereticMaterial::revertToStart(void)
{
	CrotMax = 0.0;
	CrotMin = 0.0;
	CrotPu = 0.0;
	CrotNu = 0.0;
	CenergyD = 0.0;
	CloadIndicator = 0;

	Cstress = 0.0;
	Custress = 0.0;
	Cstrain = 0.0;
	Cnlstrain = 0.0;
	CnlstrainMax = 0.0;
	CnlstrainMin = 0.0;

	Tstrain = 0;
	Tnlstrain = 0.0;
	Tstress = 0;
	Tustress = 0;
	Ttangent = E1p;
	Tutangent = E1p;

	return 0;
}

UniaxialMaterial*
HystereticMaterial::getCopy(void)
{
	HystereticMaterial *theCopy = new HystereticMaterial (this->getTag(),
		mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
		mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
		pinchX, pinchY, damfc1, damfc2, beta);

	theCopy->CrotMax = CrotMax;
	theCopy->CrotMin = CrotMin;
	theCopy->CrotPu = CrotPu;
	theCopy->CrotNu = CrotNu;
	theCopy->CenergyD = CenergyD;
	theCopy->CloadIndicator = CloadIndicator;
	theCopy->Cstress = Cstress;
	theCopy->Custress = Custress;
	theCopy->Cstrain = Cstrain;
	theCopy->Cnlstrain = Cnlstrain;
	theCopy->CnlstrainMax = CnlstrainMax;
	theCopy->CnlstrainMin = CnlstrainMin;
	theCopy->Ttangent = Ttangent;
	theCopy->Tutangent = Tutangent;

	return theCopy;
}

int
HystereticMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(32);
  
  data(0) = this->getTag();
  data(1) = mom1p;
  data(2) = rot1p;
  data(3) = mom2p;
  data(4) = rot2p;
  data(5) = mom3p;
  data(6) = rot3p;
  data(7) = mom1n;
  data(8) = rot1n;
  data(9) = mom2n;
  data(10) = rot2n;
  data(11) = mom3n;
  data(12) = rot3n;
  data(13) = pinchX;
  data(14) = pinchY;
  data(15) = damfc1;
  data(16) = damfc2;
  data(17) = beta;
  data(18) = CrotMax;
  data(19) = CrotMin;
  data(20) = CrotPu;
  data(21) = CrotNu;
  data(22) = CenergyD;
  data(23) = CloadIndicator;
  data(24) = Cstress;
  data(25) = Cstrain;
  data(26) = Ttangent;
  data(27) = Cnlstrain;
  data(28) = CnlstrainMax;
  data(29) = CnlstrainMin;
  data(30) = Custress;
  data(31) = Tutangent;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "HystereticMaterial::sendSelf() - failed to send data\n";


  return res;
}

int
HystereticMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(32);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
      opserr << "HystereticMaterial::recvSelf() - failed to receive data\n";
      return res;
  }
  else {
    this->setTag((int)data(0));
    mom1p = data(1);
    rot1p = data(2);
    mom2p = data(3);
    rot2p = data(4);
    mom3p = data(5);
    rot3p = data(6);
    mom1n = data(7);
    rot1n = data(8);
    mom2n = data(9);
    rot2n = data(10);
    mom3n = data(11);
    rot3n = data(12);
    pinchX = data(13);
    pinchY = data(14);
    damfc1 = data(15);
    damfc2 = data(16);
    beta = data(17);

    CrotMax = data(18);
    CrotMin = data(19);
    CrotPu = data(20);
    CrotNu = data(21);
    CenergyD = data(22);
    CloadIndicator = int(data(23));
    Cstress = data(24);
    Cstrain = data(25);
    Ttangent = data(26);
    Cnlstrain = data(27);
    CnlstrainMax = data(28);
    CnlstrainMin = data(29);
    Custress = data(30);
    Tutangent = data(31);

    // set the trial values
    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TrotPu = CrotPu;
    TrotNu = CrotNu;
    TenergyD = CenergyD;
    TloadIndicator = CloadIndicator;
    Tstress = Cstress;
    Tustress = Custress;
    Tstrain = Cstrain;
    Tnlstrain = Cnlstrain;
    TnlstrainMax = CnlstrainMax;
    TnlstrainMin = CnlstrainMin;
  }

  // Set envelope slopes
  this->setEnvelope();
  
  return 0;
}
    
void
HystereticMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "HHystereticMaterial, tag: " << this->getTag() << endln;
        s << "s1p: " << mom1p << endln;
        s << "e1p: " << rot1p << endln;
        s << "E1p: " << E1p << endln;
        s << "s2p: " << mom2p << endln;
        s << "e2p: " << rot2p << endln;
        s << "E2p: " << E2p << endln;
        s << "s3p: " << mom3p << endln;
        s << "e3p: " << rot3p << endln;
        s << "E3p: " << E3p << endln;
        
        s << "s1n: " << mom1n << endln;
        s << "e1n: " << rot1n << endln;
        s << "E1n: " << E1n << endln;
        s << "s2n: " << mom2n << endln;
        s << "e2n: " << rot2n << endln;
        s << "E2n: " << E2n << endln;
        s << "s3n: " << mom3n << endln;
        s << "e3n: " << rot3n << endln;
        s << "E3n: " << E3n << endln;
        
        s << "pinchX: " << pinchX << endln;
        s << "pinchY: " << pinchY << endln;
        s << "damfc1: " << damfc1 << endln;
        s << "damfc2: " << damfc2 << endln;
        s << "energyA: " << energyA << endln;
        s << "beta: " << beta << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"HystereticMaterial\", ";
        s << "\"s1p\": " << mom1p << ", ";
        s << "\"e1p\": " << rot1p << ", ";
        s << "\"E1p\": " << E1p << ", ";
        s << "\"s2p\": " << mom2p << ", ";
        s << "\"e2p\": " << rot2p << ", ";
        s << "\"E2p\": " << E2p << ", ";
        s << "\"s3p\": " << mom3p << ", ";
        s << "\"e3p\": " << rot3p << ", ";
        s << "\"E3p\": " << E3p << ", ";
        
        s << "\"s1n\": " << mom1n << ", ";
        s << "\"e1n\": " << rot1n << ", ";
        s << "\"E1n\": " << E1n << ", ";
        s << "\"s2n\": " << mom2n << ", ";
        s << "\"e2n\": " << rot2n << ", ";
        s << "\"E2n\": " << E2n << ", ";
        s << "\"s3n\": " << mom3n << ", ";
        s << "\"e3n\": " << rot3n << ", ";
        s << "\"E3n\": " << E3n << ", ";
        
        s << "\"pinchX\": " << pinchX << ", ";
        s << "\"pinchY\": " << pinchY << ", ";
        s << "\"damfc1\": " << damfc1 << ", ";
        s << "\"damfc2\": " << damfc2 << ", ";
        s << "\"energyA\": " << energyA << ", ";
        s << "\"beta\": " << beta << "}";
    }
}

void
HystereticMaterial::setEnvelope(void)
{
	E1p = mom1p/rot1p;
	E2p = (mom2p-mom1p)/(rot2p-rot1p);
	E3p = (mom3p-mom2p)/(rot3p-rot2p);

	E1n = mom1n/rot1n;
	E2n = (mom2n-mom1n)/(rot2n-rot1n);
	E3n = (mom3n-mom2n)/(rot3n-rot2n);

	Eup = E1p;
	if (E2p > Eup) Eup = E2p;
	if (E3p > Eup) Eup = E3p;

	Eun = E1n;
	if (E2n > Eun) Eun = E2n;
	if (E3n > Eun) Eun = E3n;
}

double
HystereticMaterial::posEnvlpStress(double strain)
{
  if (strain <= 0.0) {
    // threshhold the envelope to 0 in the negative strain range
    return 0.0;
  } else if (strain <= rot1p) {
    // initial elastic range
    return E1p*strain;
  } else if (strain <= rot1p + (mom2p - mom1p)/E2p) {
    // first plastic (hardening) branch
    return mom1p + E2p*(strain-rot1p);
  } else {
    // second (perfectly) plastic branch
    return mom2p;
  }
}

double
HystereticMaterial::negEnvlpStress(double strain)
{
  if (strain >= 0.0) {
    // threshhold the envelope to 0 in the positive strain range
    return 0.0;
  } else if (strain >= rot1n) {
    // initial elastic range
    return E1n*strain;
  } else if (strain >= rot1n + (mom2n - mom1n)/E2n) {
    // first plastic (hardening) branch
    return mom1n + E2n*(strain-rot1n);
  } else {
    // second (perfectly) plastic branch
    return mom2n;
  }
}

double
HystereticMaterial::posEnvlpTangent(double strain)
{
  if (strain < 0.0) {
    return E1p*1.0e-9;
  } else if (strain <= rot1p) {
    return E1p;
  } else if (strain <= rot1p + (mom2p - mom1p)/E2p) {
    return E2p;
  } else {
    return E1p*1.0e-9;
  }
}

double
HystereticMaterial::negEnvlpTangent(double strain)
{
  if (strain > 0.0) {
    return E1n*1.0e-9;
  } else if (strain >= rot1n) {
    return E1n;
  } else if (strain >= rot1n + (mom2n - mom1n)/E2n) {
    return E2n;
  } else {
    return E1n*1.0e-9;
  }
}

double
HystereticMaterial::posEnvlpRotlim(double strain)
{
  // The positive strain limit corresponds to the max allowable value of
  // the max strain value. A limiting maximum value for the max pos. strain
  // exists only when there is softening, and when the stress would otherwise go to zero.
  // Having this limit prevents the stress from softening into negative stresses.

  // set the positive strain limit to the max floating point value (+infinity)
  return POS_INF_STRAIN;
}

double
HystereticMaterial::negEnvlpRotlim(double strain)
{
  return NEG_INF_STRAIN;
}
