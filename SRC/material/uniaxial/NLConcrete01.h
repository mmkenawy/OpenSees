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
                                                                        
// $Revision: 1.0 $
// $Date: 2008-08-26 16:23:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NLConcrete01.h,v $
                                                                        
                                                                        
#ifndef NLConcrete01_h
#define NLConcrete01_h

// Written: Maha Kenawy - based on Concrete01 by MHS
// Created: August 2020
//
// Description: This file contains the class definition for 
// NLConcrete01.h adapted from Concrete01 with nonlocal modifications
//   - Modified Kent-Park envelope
//   - No tension
//   - Linear unloading/reloading


#include <UniaxialMaterial.h>

class NLConcrete01 : public UniaxialMaterial
{
 public:
  NLConcrete01 (int tag, double fpc, double eco, double fpcu, double ecu);
  NLConcrete01 ();
  ~NLConcrete01();

  const char *getClassType(void) const {return "NLConcrete01";};
  
  int setNLStrain(double nlstrain);
  int setTrialStrain(double strain, double strainRate = 0.0); 
  int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return 2.0*fpc/epsc0;}

  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int    setParameter             (const char **argv, int argc, Parameter &param);
  int    updateParameter          (int parameterID, Information &info);
  int    activateParameter        (int parameterID);
  double getStressSensitivity     (int gradIndex, bool conditional);
  int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////

  int getVariable(const char *variable, Information &);
  //by SAJalali
  double getEnergy() { return EnergyP; }

 protected:

 private:
  /*** Material Properties ***/
  double fpc;    // Compressive strength
  double epsc0;  // Strain at compressive strength
  double fpcu;   // Crushing strength
  double epscu;  // Strain at crushing strength
  
  /*** CONVERGED History Variables ***/
  double CminStrain;   // Smallest previous concrete strain (compression)
  double CunloadSlope; // Unloading (reloading) slope from CminStrain
  double CendStrain;   // Strain at the end of unloading from CminStrain
  
  /*** CONVERGED State Variables ***/
  double Cstrain;
  double Cstress;   
  double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
  // Storing it is better than recomputing it!!!
  double Cnlstrain;
  
  /*** TRIAL History Variables ***/
  double TminStrain;
  double TunloadSlope;
  double TendStrain;
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
  double Tnlstrain;
  
  void determineTrialState (double dStrain);
  
  void reload();
  void unload();
  void envelope();
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Matrix *SHVs;
  // AddingSensitivity:END ///////////////////////////////////////////

  //by SAJalali
  double EnergyP;
};


#endif


