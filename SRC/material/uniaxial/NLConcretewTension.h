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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NLConcretewTension.h,v $

#ifndef NLConcretewTension_h
#define NLConcretewTension_h

// Written: Maha Kenawy
// Created: August 2020
//
// refer to paper:
//

#include <UniaxialMaterial.h>
#include <Matrix.h>


class NLConcretewTension : public UniaxialMaterial
{
  public:
    NLConcretewTension(int tag, double E, double fc,
		      double ec0, double Ed, double ft, double eft);
    ~NLConcretewTension();

    const char *getClassType(void) const {return "NLConcretewTension";};

    int setTrialStrain(double strain, double strainRate = 0.0);
    int setNLStrain(double nlstrain);
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

  protected:
    
  private:
    // Material parameters
    double E;	// Elastic modulus
    double fc;	// compressive strength
    double ec0;	// yield strain in compression
    double Ed;	// softening slope
    double ft; // tensile strength
    double eft; // tensile strain at zero tensile stress
	
    // Committed history variables
    double CPlasticStrain;	// Committed plastic strain
    double CDam;
    double Ckk;
    double Ckc;		// Committed accumulated plastic strain
    double Ckt;

	// Trial history variables
    double TPlasticStrain;	// Trial plastic strain
    double TDam;
    double Tkk;
    double Tkc;		// Trial accumulated plastic strain
    double Tkt;

    // Trial state variables
    double TStrain;		// Trial strain
    double Tnlstrain;
    double TStress;		// Trial stress
    double Tefstress;
    double TTangent;	// Trial tangent
    double TDc;
    double TDt;
    double Tkdc1;
    double Tkdt1;

    double CStrain;		// Committed strain
    double Cnlstrain;
    double CStress;		// Committed stress
    double Cefstress;
    double CTangent;	// Committed tangent
    double CDc;
    double CDt;
    double Ckdc1;
    double Ckdt1;
};


#endif
