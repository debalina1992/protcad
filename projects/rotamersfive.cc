#include <iostream>
#include <string>
#include <time.h>
#include <sstream>
#include <stdio.h>
#include "ensemble.h"
#include "PDBInterface.h"
#include "assert.h"
#include <string.h>
#include <vector>
#include "typedef.h"
#include "ran.h"
#include "molecule.h"
#include "chain.h"
#include "chainModBuffer.h"

int main (int argc, char* argv[])
{
if (argc !=2)
	{
		cout << "rotamers <input.pdb>" << endl;
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
	//string outFile;
	
UInt randchain, randres,size1,size2,size3,size4,size5,randrestype;
UInt resNum, chainNum = _prot->getNumChains();
UIntVec allowedRots;
UIntVec allowedRots2;
UIntVec allowedRots3,allowedRots4,allowedRots5;

cout << "INDEX A.B  == A: Mainchain Brainch Point Conf :: B: Sidechain Branch Point Conf" << endl;
		
		int a=0;

		randchain = rand() % chainNum;
		resNum = _prot->getNumResidues(randchain);
		randres = rand() % resNum;
		randrestype = _prot->getTypeFromResNum(randchain, randres);

allowedRots = _prot->getAllowedRotamers(randchain, randres, randrestype, 0);
size1 = allowedRots.size();
for (UInt i = 0; i < size1; i ++)
	{ 
	_prot->setRotamerWBC(randchain, randres, 0, allowedRots[i]);
	allowedRots2 = _prot->getAllowedRotamers(randchain, randres, randrestype, 1);
 	size2 = allowedRots2.size();
		for(UInt j = 0; j < size2; j ++)
		{	
		_prot->setRotamerWBC(randchain, randres, 1, allowedRots2[j]);
		allowedRots3 = _prot->getAllowedRotamers(randchain, randres, randrestype, 2);
 		size3 = allowedRots3.size();
		for(UInt k = 0; k < size3; k ++)
			{	
			_prot->setRotamerWBC(randchain, randres, 2, allowedRots3[k]);
			allowedRots4 = _prot->getAllowedRotamers(randchain, randres, randrestype, 3);
 			size4 = allowedRots4.size();
			for(UInt l = 0; l < size4; l ++)
				{	
				_prot->setRotamerWBC(randchain, randres, 3, allowedRots4[l]);
				allowedRots5 = _prot->getAllowedRotamers(randchain, randres, randrestype, 4);
 				size5 = allowedRots5.size();
				for(UInt m = 0; m < size5; m ++)
					{	
					_prot->setRotamerWBC(randchain, randres, 4, allowedRots5[m]);

					a++;				
					ostringstream convert; 
					convert << a;    
					string Result = convert.str();
					string outFile= "conf"+Result+".pdb";

					pdbWriter(_prot, outFile);
					double intra = _prot->intraSoluteEnergy(true);
					cout << "The energy for conformation " << a << " : " << i << "." << j << "." << k << "." << l << "." << m << " is " << intra << " FILE: " << outFile << endl;	
					}
				}
			}
		}
	}
return 0;
}
