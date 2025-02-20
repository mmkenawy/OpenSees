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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
                                                                        
#ifndef NLTruss_h
#define NLTruss_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for NLTruss. A NLTruss object
// provides the abstraction of the small deformation bar element. Each NLtruss
// object is associated with a material object. This NLTruss element will work
// in 1d, 2d or 3d problems.
//
// What: "@(#) NLTruss.h, revA"

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class NLTruss : public Element
{
  public:
    NLTruss(int tag, int dimension,
	  int Nd1, int Nd2, 
	  UniaxialMaterial &theMaterial,
	  double A, double nllength = 0.0, double rho = 0.0,
	  int doRayleighDamping = 0,
      int cMass = 0);
    
    NLTruss();
    ~NLTruss();

    const char *getClassType(void) const {return "NLTruss";};

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    int computeNLStrain(void);
    
    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getKi(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);    
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int		   addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double computeCurrentStrain(void) const;
    double computeCurrentStrainRate(void) const;
    
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;  // pointer to a material
    ID  connectedExternalNodes;     // contains the tags of the end nodes
    int dimension;                  // NLtruss in 2 or 3d domain
    int numDOF;	                    // number of dof for NLtruss

    Vector *theLoad;    // pointer to the load vector P
    Matrix *theMatrix;  // pointer to objects matrix (a class wide Matrix)
    Vector *theVector;  // pointer to objects vector (a class wide Vector)

    double L;               // length of NLtruss based on undeformed configuration
    double A;               // area of NLtruss
    double nllength;        // nllength parameter
    double rho;             // rho: mass density per unit length
    int doRayleighDamping;  // flag to include Rayleigh damping
    int cMass;              // consistent mass flag

    double cosX[3];  // direction cosines

    Node *theNodes[2];
    double *initialDisp;

	
// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    Vector *theLoadSens;
// AddingSensitivity:END ///////////////////////////////////////////

    // static data - single copy for all objects of the class	
    static Matrix NLtrussM2;   // class wide matrix for 2*2
    static Matrix NLtrussM4;   // class wide matrix for 4*4
    static Matrix NLtrussM6;   // class wide matrix for 6*6
    static Matrix NLtrussM12;  // class wide matrix for 12*12
    static Vector NLtrussV2;   // class wide Vector for size 2
    static Vector NLtrussV4;   // class wide Vector for size 4
    static Vector NLtrussV6;   // class wide Vector for size 6
    static Vector NLtrussV12;  // class wide Vector for size 12
};

#endif
