#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "Node.h"
#include "Group.h"
#include "Element.h"
#include "Material.h"
#include "Prescribe.h"
#include "Load.h"
#include "mkl.h"
#include "gidpost.h"
#include "mpi.h"

using namespace std;
class FEMOOP
{
	ofstream chk, res;
	time_t  tNow;
	struct tm   tmLocal;
	string probn, text;
	int  nnode, ndim, nelem, ngroup, nmat, nstep;
	string GidResFile;
	bool Rebuild;
	int myid, numprocs;

	double *A, Error,ErrorSum, *Result;

	int TotalDegreeOfFreedom, *DegreeOfFreedom;
	double ** GlobalStiffMatrix, **LoadMatrix, **RightHand, *InitialForce;
	double **InitialDisplace, *DispLoad,*InteractLoad;
public:
	int InputFile();
	int SendID(int numprocs, int myid);
	int SetProbn(string probn);
	int OpenCheckFile();
	int ShowPassTime();
	int ReadControl();
	int ReadAllFile();
	int GIDOutMesh();
	int TotalDOF(int istep,bool Rebuild);
	int GlobalStiff();
	int ApplyLoad(int istep);
	int Solve(int istep);
	int StrainForce();
	int GIDOutResult(int istep);
	int StepCycle();
	//

	//! PostProcess().
	/*!
	
	*/
	int PostProcess(int istep);
};


