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
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);
	//string outFile;

UInt size2,size3,size4;
UIntVec allowedRots2;
UIntVec allowedRots3,allowedRots4,allowedRots5;
ofstream fout;


fout.open("energy.dat",ios::app);    // open file for appending
 
UInt a=0;
UInt angle=0;
allowedRots2 = _prot->getAllowedRotamers(0, 31,54, 1);
size2 = allowedRots2.size();
allowedRots3 = _prot->getAllowedRotamers(0, 31, 54, 2);
 size3 = allowedRots3.size();
allowedRots4 = _prot->getAllowedRotamers(0, 31, 54, 3);
 size4 = allowedRots4.size();

	for (UInt i = 0; i < 36; i ++)
	{  		_prot->setChi(0, 31, 0, 0, angle);
			_prot->setChi(1, 31, 0, 0, angle);
			_prot->setChi(2, 31, 0, 0, angle);
			_prot->setChi(3, 31, 0, 0, angle);
			angle = angle+10;
		for(UInt j = 0; j < size2; j ++)
		{			
			for(UInt k = 0; k < size3; k ++)
			{				
				for(UInt l = 0; l < size4; l ++)
				{	
					
						#pragma omp parallel for 
						for(UInt p = 0; p < 4 ; p ++)
						{		
						_prot->setRotamerWBC(p, 31, 1, allowedRots2[j]);
	  					_prot->setRotamerWBC(p, 31, 2, allowedRots3[k]);
						_prot->setRotamerWBC(p, 31, 3, allowedRots4[l]);
						}
						
						a++;				
						ostringstream convert; 
						convert << a;    
						string Result = convert.str();
						string outFile= "conf"+Result+".pdb";
	
						pdbWriter(_prot, outFile);
						double intra = _prot->intraSoluteEnergy(false);
						fout << "chi=" << "	"<< angle-10 << "	" << j << "." << k << "." << l << "	" << intra << "  " << outFile << endl;		

					//}
				}
			}
		}
	}
return 0;
}
