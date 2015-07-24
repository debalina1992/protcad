#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include <stdio.h>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[])
{
// if (argc !=7)
//	{
//		cout << "mutation4 <input.pdb> <output1.pdb> <output2.pdb> <output3.pdb> <output4.pdb> <output5.pdb>" << endl;
//		exit(1);
//	}
	enum aminoAcid {A,R,N,D,Dh,C,Cx,Q,E,Eh,Hd,He,Hn,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dQ,dE,dEh,dHd,dHe,dHn,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,HCE,PCH};
	string infile=argv[1];
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* _prot = static_cast<protein*>(pMol);
	_prot->silenceMessages();
	residue::setCutoffDistance(9.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
    	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(1.0);
	solvation::setItsScaleFactor(0.0);
	string outFile;

	 	//UInt resNum = _prot->getNumResidues(0);
		//UInt randres = rand() % resNum;
		//_prot->activateForRepacking(0,0);
		//_prot->mutateWBC(0,0,19);
		pdbWriter(_prot, "mutatedoutput.pdb");
return 0;
}

