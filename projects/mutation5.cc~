#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include <stdio.h>
#include "ensemble.h"
#include "PDBInterface.h"
int main (int argc, char* argv[2])
{ if (argc !=6)
	{
		cout << "mutation4 <input.pdb> <output1.pdb> <output2.pdb> <output3.pdb> <output4.pdb> <output5.pdb>" << endl;
		exit(1);
	}
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


/*string residue=argv[2];

//string to char conversion:


char *res = new char[residue.length() + 1];
strcpy(res, residue.c_str());


// do stuff
//delete [] cstr;


//	*/

int a=1;
{ while ( a<6) {
		double intra = _prot->intraSoluteEnergy(true);
		
	 	UInt resNum = _prot->getNumResidues(0);
		UInt randres = rand() % resNum;
		UInt restype = _prot->getTypeFromResNum(0,randres);
		

		_prot->mutateWBC(0,randres,R);
		_prot->protOptSolvent(500);
		double intra2=_prot->intraSoluteEnergy(true);
	if (intra2 < intra)
		{outFile = argv[a+1];
		a++;		
		pdbWriter(_prot, outFile);
		} 
else 		{_prot->mutate(0,randres,restype);
		//undoLastMutation()
		}
	}
}
return 0;
}

