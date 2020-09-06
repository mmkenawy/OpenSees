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

// Written: QGU. UCSD Group for response sensitivity analysis
//  Contact: Quan Gu(guquan2005@gmail.com)  
//           Joel P. Conte (jpconte@ucsd.edu)
//           Michele Barbato (mbarbato@lsu.edu)
//
// Created: May 2009

// refer to paper:
// Conte, J. P., Vijalapura, P., and Meghella, M., "Consistent Finite Element Response Sensitivities Analysis," 
// Journal of Engineering Mechanics, ASCE, Vol. 129, No. 12, pp. 1380-1393, 2003.
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



NLConcrete::NLConcrete(int pTag, double pE, double pSigmaY,
				     double pHkin, double pHiso)
:UniaxialMaterial(pTag,MAT_TAG_NLConcrete),
 E(pE), fc(pSigmaY), ec0(pHiso), Ed(pHkin)
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	// Initialize variables
    this->revertToStart();
}

NLConcrete::~NLConcrete()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
	SHVs =0;
// AddingSensitivity:END //////////////////////////////////////
}

int NLConcrete::setNLStrain(double nlstrain)
{
    Tnlstrain = nlstrain;

    return 0;
}

int 
NLConcrete::setTrialStrain (double strain, double strainRate)
{

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

        // for continued plasticity model (will not be properly calibrated)
       // double kf = (E*ef - fy)/(H+E);

        double alphalin = 1.0;

//        double mu = 0.85;
//        double H2 = (E*ef*(1 - mu) - fc)/(mu*ef - k0);
//        if (H2 < 0.0)
//                H2 = 0.01;
        double H2 = 0.15*E;


//      // damage model parameters
//        double e1 = fy/E;
//            double logterm1 = log10(fc/(E*e1));
//            double logterm2 = log10((ef-ec0)/(ef-e1)) - (ec0/(ef - ec0))*log10(e1/ec0);
//            double alpha = logterm1/logterm2;
//            double beta = 1 - (alpha*ec0/(ef - ec0));


    // ------ Elastic trial -------


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

                //////////// beginning of plasticity part - comment out for damage-only model ////////////////

        double TrialStresstmp = E*(TStrain - CPlasticStrain);
        double TStresstmp = TrialStresstmp;

         if (TrialStresstmp > 0.0)
    {
        TStress = 0.0;
        TTangent = 1.0e-10;
//        TTangent = E;
    }

        else {

             if (TAccumulatedPlasticStrain <= k0)
             {
                Mod = H;
                Tystress = fy + Mod*TAccumulatedPlasticStrain;
             }
             else {
                Mod = H2;
                Tystress = fc + Mod*(TAccumulatedPlasticStrain - k0);
             }

    // Compute trial stress relative to committed back stress

    // Compute yield criterion
    double f = fabs(TrialStresstmp) - Tystress;

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
		//if (f <= 1.0e-8) {
		// Set trial tangent
		TTangent = (1 - TD)*E;
                elastic = 1;
    }

    // ------- Plastic step ... perform return mapping algorithm ---
    else {
        elastic = 0;
      deltaLambda = f/(E+Mod);
//        // alternative method for finding plastic strain
//      int numiter = 3;
//      Matrix x(3,numiter+1);
//      Matrix R(3,numiter+1);
//      Matrix Jinv(3,3);
//      Vector dx(3);
//      R.Zero();
//      x.Zero();
////      //R(3,0) = TStresstmp - Tystress;
//      int sign = (TrialStresstmp < 0) ? -1 : 1;
//      for (int i = 0; i < numiter; i++) {
//          R(0,i) = x(0,i) + E*x(2,i)*sign - TrialStresstmp;
//          R(1,i) = x(1,i) - CAccumulatedPlasticStrain - x(2,i);
//          if (x(1,i) <= k0)
//            //R(2,i) = fabs(x(0,i)) - (fy + H*x(1,i));
//              R(2,i) = fabs(x(0,i)) - (fc/3)*(1 + 4*x(1,i)/k0 - 2*pow(x(1,i),2)/pow(k0,2));
//             else
//                 R(2,i) = fabs(x(0,i)) - fc;
//
//          // calculate the inverse of the jacobian
//          double C = 1/(3*E*pow(k0,2) + 4*fc*k0 - 4*fc*x(1,i));
//          //double C = 1/(E+H);
//
//          Jinv(0,0) = C*-4*fc*(x(1,i) - k0);
//          //Jinv(0,0) = C*H;
//          Jinv(0,1) = -Jinv(0,0)*E;
//          Jinv(0,2) = C*-3*E*pow(k0,2);
//  //        Jinv(0,2) = C*-E;
//          Jinv(1,0) = C*-3*pow(k0,2);
//          //Jinv(1,0) = -C;
//          Jinv(1,1) = -Jinv(0,2);
//
//        if (x(1,i) > k0) {
//          Jinv(0,0) = 0.0;
//          Jinv(0,1) = 0.0;
//          Jinv(0,2) = -1;
//          Jinv(1,0) = -1/E;
//          Jinv(1,1) = 1;
//        }
//          Jinv(1,2) = Jinv(1,0);
//          Jinv(2,0) = Jinv(1,0);
//          Jinv(2,1) = -Jinv(0,0);
//          Jinv(2,2) = Jinv(1,0);
////
////          dx(0) = -1/(E+H)*(H*R(0,i) - E*H*R(1,i) - E*R(2,i));
////          dx(1) = -1/(E+H)*(-1*R(0,i) + E*R(1,i) - 1*R(2,i));
////          dx(2) = -1/(E+H)*(-1*R(0,i) - H*R(1,i) - 1*R(2,i));
//          dx(0) = -Jinv(0,0)*R(0,i) - Jinv(0,1)*R(1,i) - Jinv(0,2)*R(2,i);
//          dx(1) = -Jinv(1,0)*R(0,i) - Jinv(1,1)*R(1,i) - Jinv(1,2)*R(2,i);
//          dx(2) = -Jinv(2,0)*R(0,i) - Jinv(2,1)*R(1,i) - Jinv(2,2)*R(2,i);
//          x(0,i+1) = x(0,i) + dx(0);
//          x(1,i+1) = x(1,i) + dx(1);
//          x(2,i+1) = x(2,i) + dx(2);
////
//      }
//      TStresstmp = x(0,numiter);
//      deltaLambda = x(2,numiter);
//      TAccumulatedPlasticStrain = x(1,numiter);

      // Find sign of xsi
      int sign = (TStresstmp < 0) ? -1 : 1;

	  TPlasticStrain = CPlasticStrain +deltaLambda*sign;

	  TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;

	  TStresstmp = TrialStresstmp - E*deltaLambda*sign;
          //TTangent = E*H/(E + H);
          double estraininc = dStrain - deltaLambda*sign;
	
// calculate the damage
          //////// new damage model
                   // calculate the nonlocal damage variables kd1 and kd2
          double dnlk = std::max(fabs(dnlstrain) - fabs(estraininc),0.0);
          //double dnlk = std::max(deltaLambda + m*fabs(dnlstrain) - m*fabs(dStrain),0.0);
          Tnlk = Cnlk + dnlk;
          //Tnlk = TAccumulatedPlasticStrain + m*fabs(Tnlstrain) - m*fabs(TStrain);
          if (Tnlk > k0)
          {
              double kd1 = Tnlk - k0;
//              Tkd2 = fabs(Tstresstmp/E);
//              Tkd2 = (fc + Mod*kd1)/E;
//              if (Tkd2 < Ckd2)
//                  Tkd2 = Ckd2;
              //TD = (E*efs*Tkd2 + fc*kd1 - fc*efs)/(Tkd2*(E*efs - fc));
              TD = (E*kd1*(efs*H2 + fc))/((fc + H2*kd1)*(E*efs - fc));
              dDoverdk = fc*E*(efs*H2 + fc)/(E*efs - fc)/pow((fc + H2*kd1),2);
          }

          if (TD > 1.0)
                  TD = 1.0;
          if (fabs(TD - CD) < 1.0e-8)
              dDoverdk = 0.0;


          ////// old damage model
//        double dnlk  = std::max(deltaLambda + m*fabs(dnlstrain) - m*fabs(dStrain),0.0);
//          Tnlk = Cnlk + dnlk;
//          //Tnlk = TAccumulatedPlasticStrain + m*fabs(Tnlstrain) - m*fabs(TStrain);
//          if (Tnlk > kf)
//              TD = 1.0;
//          else if (Tnlk > k0) {
//          //TD = std::max(1 - exp(-(Tnlk-k0)/kf),0.0);
//          //TD = 0.0;
//          TD = std::min(std::max(1 - pow((kf - Tnlk)/(kf - k0),alphalin),0.0),1.0);
//
//          if (TD > 0.0)
//            //dDoverdk = 1/(kf-k0);
//              dDoverdk = alphalin*pow((kf - Tnlk)/(kf - k0),alphalin)/(kf - Tnlk);
//          //dDoverdk = exp(-(Tnlk-k0)/kf)/kf;
//          }
//          if (TD > 0.0)
//                  TTangent = -fc/efs;
         //TTangent = ((1 - TD)*Mod/(E+Mod) - TStresstmp*dDoverdk*-1/(E+Mod))*E;
        //TTangent = E;
        TTangent = ((1 - TD)*Mod/(E+Mod) - TStresstmp*dDoverdk*-1/(E+Mod))*E;
        if (TD > 0.999)
            TTangent = 1.0e-8;
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

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

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
    s << "  sigmaY: " << fc << endln;
    s << "  Hiso: " << ec0 << endln;
    s << "  Hkin: " << Ed << endln;
   
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
NLConcrete::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(2, this);

  if ((strcmp(argv[0],"H_kin") == 0)||(strcmp(argv[0],"Hkin") == 0))
    return param.addObject(3, this);

  if ((strcmp(argv[0],"H_iso") == 0)||(strcmp(argv[0],"Hiso") == 0))
    return param.addObject(4, this);

  return -1;
}

int
NLConcrete::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->fc = info.theDouble;
		break;
	case 2:
		this->E = info.theDouble;
		break;
	case 3:
		this->ec0 = info.theDouble;
		break;
	case 4:
		this->Ed = info.theDouble;
		break;
	default:
		return -1;
	}

	return 0;
}



int
NLConcrete::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}




double
NLConcrete::getStrainSensitivity(int gradIndex)
{
//
//	if (SHVs ==0) return 0.0;
//	else{
//		double sensitivity =(*SHVs)(4,gradIndex-1); // unconditional stress sensitivity
		return 0.0;
	}
//}

double
NLConcrete::getStressSensitivity(int gradIndex, bool conditional)
{
//
//	if (conditional == false) {  // return stress sensitivity for recorder purpose
//		if (SHVs ==0) return 0.0;
//		else {
//			double sensitivity =(*SHVs)(3,gradIndex-1); // unconditional stress sensitivity
//			return sensitivity;
//		}
//	}
//
//	// First set values depending on what is random
//	double SigmaYSensitivity = 0.0;
//	double ESensitivity = 0.0;
//	double HkinSensitivity = 0.0;
//	double HisoSensitivity = 0.0;
//
//	if (parameterID == 1) {  // sigmaY
//		SigmaYSensitivity = 1.0;
//	}
//	else if (parameterID == 2) {  // E
//		ESensitivity = 1.0;
//	}
//	else if (parameterID == 3) {  // Hkin
//		HkinSensitivity = 1.0;
//	}
//	else if (parameterID == 4) {  // Hiso
//		HisoSensitivity = 1.0;
//	}
//	else {
//		// Nothing random here, but may have to return something in any case
//	}
//
//	double TStrainSensitivity = 0.0;
//
//	// Then pick up history variables for this gradient number
//	double CPlasticStrainSensitivity = 0.0;
//	double CBackStressSensitivity	 = 0.0;
//	double CAccumulatedPlasticStrainSensitivity	 = 0.0;
//	if (SHVs != 0) {
//		CPlasticStrainSensitivity = (*SHVs)(0,gradIndex);
//		CBackStressSensitivity	 = (*SHVs)(1,gradIndex);
//		CAccumulatedPlasticStrainSensitivity = (*SHVs)(2,gradIndex);
//	}
//
//    // ------ Elastic trial -------
//
//	double deltaLambda = 0.0;
//	TPlasticStrain = CPlasticStrain;
//	TBackStress = CBackStress;
//	TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
//    TStress = E * (TStrain-CPlasticStrain);
//	double TStressSensitivity = E*(TStrainSensitivity-CPlasticStrainSensitivity)+ESensitivity*(TStrain-CPlasticStrain);
//	double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;
//
//    // Compute trial stress relative to committed back stress
//    double xsi = TStress - TBackStress;
//
//    // Compute yield criterion
//    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);
//
//	double sensitivity;
//
//
//    // Elastic step ... no updates required
//    if (f <= -DBL_EPSILON * E) {
//		//if (f <= 1.0e-8) {
//		// Set trial tangent
//		TTangent = E;
//		sensitivity = TStressSensitivity;
//	}
//
//    // ------- Plastic corrector ... perform return mapping algorithm ---
//    else {
//      deltaLambda = (fabs(xsi)-CSigmaY)/(E+Hkin+Hiso);
//
//	  // Find sign of xsi
//      int sign = (xsi < 0) ? -1 : 1;
//	  TPlasticStrain = CPlasticStrain +deltaLambda*sign;
//
//	  TBackStress = CBackStress + Hkin*deltaLambda*sign;
//
//	  TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;
//
//	  TStress = E * (TStrain - TPlasticStrain);
//
//      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);
//
//	  double CSigmaYSensitivity = SigmaYSensitivity + HisoSensitivity*CAccumulatedPlasticStrain
//		   + Hiso*CAccumulatedPlasticStrainSensitivity;
//
//	  double deltaLambdaSensitivity = ((TStressSensitivity-CBackStressSensitivity)*sign-CSigmaYSensitivity)/(E+Hkin+Hiso)
//		  - (ESensitivity+HkinSensitivity+HisoSensitivity)*((E * (TStrain-CPlasticStrain)-CBackStress)*sign-CSigmaY)/pow((E+Hkin+Hiso),2.0);
//
//	  double TPlasticStrainSensitivity = CPlasticStrainSensitivity + deltaLambdaSensitivity*sign;
//
//	  sensitivity = E * (TStrainSensitivity - TPlasticStrainSensitivity)+ESensitivity * (TStrain - TPlasticStrain);
//
//   }
	return 0.0;
}
//
//
//
double
NLConcrete::getInitialTangentSensitivity(int gradIndex)
{
//	// For now, assume that this is only called for initial stiffness
//	if (parameterID == 2) {
//		return 1.0;
//	}
//	else {
		return 0.0;
	}
//}
//
//
int
NLConcrete::commitSensitivity(double TStrainSensitivity, int gradIndex, int numGrads)
{
//
//
//	if (SHVs == 0) {
//		SHVs = new Matrix(5,numGrads);
//		SHVs->Zero();
//	}
//
//
//
//	// First set values depending on what is random
//	double SigmaYSensitivity = 0.0;
//	double ESensitivity = 0.0;
//	double HkinSensitivity = 0.0;
//	double HisoSensitivity = 0.0;
//
//	if (parameterID == 1) {  // sigmaY
//		SigmaYSensitivity = 1.0;
//	}
//	else if (parameterID == 2) {  // E
//		ESensitivity = 1.0;
//	}
//	else if (parameterID == 3) {  // Hkin
//		HkinSensitivity = 1.0;
//	}
//	else if (parameterID == 4) {  // Hiso
//		HisoSensitivity = 1.0;
//	}
//	else {
//		// Nothing random here, but may have to return something in any case
//	}
//
//
//	// Then pick up history variables for this gradient number
//
//	double CPlasticStrainSensitivity = (*SHVs)(0,gradIndex);
//	double CBackStressSensitivity	 = (*SHVs)(1,gradIndex);
//	double CAccumulatedPlasticStrainSensitivity = (*SHVs)(2,gradIndex);
//
//
//    // ------ Elastic trial -------
//
//	double deltaLambda = 0.0;
//	TPlasticStrain = CPlasticStrain;
//	TBackStress = CBackStress;
//	TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
//    TStress = E * (TStrain-CPlasticStrain);
//	double TStressSensitivity = E*(TStrainSensitivity-CPlasticStrainSensitivity)+ESensitivity*(TStrain-CPlasticStrain);
//	double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;
//
//    // Compute trial stress relative to committed back stress
//    double xsi = TStress - TBackStress;
//
//    // Compute yield criterion
//    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);
//
//	double sensitivity;
//
//
//    // Elastic step ... no updates required
//    if (f <= -DBL_EPSILON * E) {
//		//if (f <= 1.0e-8) {
//		// Set trial tangent
//		TTangent = E;
//		sensitivity = TStressSensitivity;
//	}
//
//    // ------- Plastic step ... perform return mapping algorithm ---
//    else {
//      deltaLambda = (fabs(xsi)-CSigmaY)/(E+Hkin+Hiso);
//
//	  // Find sign of xsi
//      int sign = (xsi < 0) ? -1 : 1;
//	  TPlasticStrain = CPlasticStrain +deltaLambda*sign;
//
//	  TBackStress = CBackStress + Hkin*deltaLambda*sign;
//
//	  TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;
//
//	  TStress = E * (TStrain - TPlasticStrain);
//
//      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);
//
//	  double CSigmaYSensitivity = SigmaYSensitivity + HisoSensitivity*CAccumulatedPlasticStrain
//		   + Hiso*CAccumulatedPlasticStrainSensitivity;
//
//	  double deltaLambdaSensitivity = ((TStressSensitivity-CBackStressSensitivity)*sign-CSigmaYSensitivity)/(E+Hkin+Hiso)
//		  - (ESensitivity+HkinSensitivity+HisoSensitivity)*((E * (TStrain-CPlasticStrain)-CBackStress)*sign-CSigmaY)/pow((E+Hkin+Hiso),2.0);
//
//	  double TPlasticStrainSensitivity = CPlasticStrainSensitivity + deltaLambdaSensitivity*sign;
//
//	  sensitivity = E * (TStrainSensitivity - TPlasticStrainSensitivity)+ESensitivity * (TStrain - TPlasticStrain);
//
//	  double TAccumulatedPlasticStrainSensitivity = CAccumulatedPlasticStrainSensitivity + deltaLambdaSensitivity;
//
//	  double TBackStressSensitivity = CBackStressSensitivity + HkinSensitivity*deltaLambda*sign + Hkin*deltaLambdaSensitivity*sign;
//
//
//	  (*SHVs)(0,gradIndex) = TPlasticStrainSensitivity;
//	  (*SHVs)(1,gradIndex) = TBackStressSensitivity;
//	  (*SHVs)(2,gradIndex) = TAccumulatedPlasticStrainSensitivity;
//	  (*SHVs)(3,gradIndex) = sensitivity;      // for recorder purpose
//	  (*SHVs)(4,gradIndex) = TStrainSensitivity;  // for recorder purpose
//
//
//
//
//    }
//
//
//
 	return 0;
}
//// AddingSensitivity:END /////////////////////////////////////////////
