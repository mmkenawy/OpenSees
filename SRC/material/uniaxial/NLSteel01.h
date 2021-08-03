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
                                                                        
// $Revision: 1.14 $
// $Date: 2008-08-26 16:35:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NLSteel01.h,v $
                                                                        
                                                                        
#ifndef NLSteel01_h
#define NLSteel01_h

// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class definition for 
// NLSteel01.h
// 
//
//
// What: "@(#) NLSteel01.h, revA"


#include <UniaxialMaterial.h>

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define STEEL_01_DEFAULT_A1        0.0
#define STEEL_01_DEFAULT_A2       55.0
#define STEEL_01_DEFAULT_A3        0.0
#define STEEL_01_DEFAULT_A4       55.0

class NLSteel01 : public UniaxialMaterial
{
  public:
    NLSteel01(int tag, double fy, double E0, double b,
	   double initpStrain = 0.0, double initnStrain = 0.0, double maxpStrain = 0.0,
	   double maxnStrain = 0.0, double maxpDamage = 1.0, double maxnDamage = 1.0,
	   double a1 = STEEL_01_DEFAULT_A1, double a2 = STEEL_01_DEFAULT_A2,
	   double a3 = STEEL_01_DEFAULT_A3, double a4 = STEEL_01_DEFAULT_A4);
    NLSteel01();
    ~NLSteel01();

    const char *getClassType(void) const {return "NLSteel01";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setNLStrain(double nlstrain);
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E0;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
//// AddingSensitivity:BEGIN //////////////////////////////////////////
//    int setParameter(const char **argv, int argc, Parameter &param);
//    int    updateParameter          (int parameterID, Information &info);
//    int    activateParameter        (int parameterID);
//    double getStressSensitivity     (int gradIndex, bool conditional);
//    double getInitialTangentSensitivity(int gradIndex);
//    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
//    // AddingSensitivity:END ///////////////////////////////////////////
	//by SAJalali
	virtual double getEnergy() { return Energy; }

 protected:
    
 private:
	 double Energy;	//by SAJalali
	/*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double a1;
    double a2;
    double a3;
    double a4;  // a1 through a4 are coefficients for isotropic hardening
    double initpStrain;
    double initnStrain;
    double maxpStrain;
    double maxnStrain;
    double maxpDamage;
    double maxnDamage;
    
    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cnlstrain;
    double Custress;
    double Cstress;
    double Cpdamage;
    double Cndamage;
    double Cutangent;
    double Ctangent;

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tnlstrain;
    double Tustress;
    double Tstress;
    double Tpdamage;
    double Tndamage;
    double Tutangent;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);

    // Apply damage to the current stress state
    void applyDamage (void);

    // Determines if a load reversal has occurred based on the trial strain
    void detectLoadReversal (double dStrain);

//// AddingSensitivity:BEGIN //////////////////////////////////////////
//    int parameterID;
//	Matrix *SHVs;
//// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
