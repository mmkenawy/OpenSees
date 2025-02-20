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

// Default values for initial and max strain (used by non-local damage model) -bdg
#define NLSTEEL_01_DEFAULT_INIT_STRAIN 0.01
#define NLSTEEL_01_DEFAULT_MAX_STRAIN  0.1

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define NLSTEEL_01_DEFAULT_A1        0.0
#define NLSTEEL_01_DEFAULT_A2       55.0
#define NLSTEEL_01_DEFAULT_A3        0.0
#define NLSTEEL_01_DEFAULT_A4       55.0

class NLSteel01 : public UniaxialMaterial
{
  public:
    NLSteel01(int tag, double fy, double E0, double b,
       double initStrain = NLSTEEL_01_DEFAULT_INIT_STRAIN, double maxStrain = NLSTEEL_01_DEFAULT_MAX_STRAIN,
       double a1 = NLSTEEL_01_DEFAULT_A1, double a2 = NLSTEEL_01_DEFAULT_A2,
       double a3 = NLSTEEL_01_DEFAULT_A3, double a4 = NLSTEEL_01_DEFAULT_A4);
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
    
// AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////
	//by SAJalali
	virtual double getEnergy() { return Energy; }

 protected:
    
 private:
	 double Energy;	//by SAJalali
	/*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double initStrain; // initial strain at damage in compression for non-local buckling model -bdg
    double maxStrain;  // maximum strain at damage in compression for non-local buckling model -bdg
    double a1;
    double a2;
    double a3;
    double a4;  // a1 through a4 are coefficients for isotropic hardening

    /*** CONVERGED State Variables ***/
    double Cstrain;   // strain
    double Cnlstrain; // non-local strain
    double Cstress;   // stress
    double Ctangent;  // material stiffness (tangent)
    double Cpstrain;  // total plastic strain
    double Ceps;      // equivalent plastic strain
    double Cdamage;   // damage variable controlling compressive softening
    
    /*** TRIAL State Variables ***/
    double Tstrain;   // strain
    double Tnlstrain; // non-local strain
    double Tstress;   // stress
    double Ttangent;  // material stiffness (tangent)
    double Tpstrain;  // total plastic strain
    double Teps;      // equivalent plastic strain
    double Tdamage;   // damage variable controlling compressive softening

    // Calculates the trial state variables based on the trial strain
    void determineTrialState ();

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
