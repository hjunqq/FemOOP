#include "FEMOOP.h"
extern Node Nodes;
extern Element Elems;
extern Group Groups;
extern Load *Loads;
extern Material Mats;
extern Prescribe *Pres;
//int FEMOOP::InputFile()
//{
//	ifstream inp("1.inp");
//	if (inp.fail())
//	{
//		cout << "The input file does not exist" << endl;
//		cout << "Exit program" << endl; 
//		return -1;
//	}
//	getline(inp, text);
//	getline(inp, probn);
//	inp.close();
//	return 0;
//}
int FEMOOP::SendID(int numprocs, int myid)
{
	this->numprocs = numprocs;
	this->myid = myid;
	return 0;
}
int FEMOOP::SetProbn(string probn)
{
	this->probn = probn;
	cout << "Process " << myid << endl;
	return 0;
}
int FEMOOP::ReadControl()
{
	ifstream glb(probn + ".glb");
	if (glb.fail())
	{
		cout << "The input file does not exist" << endl;
		cout << "Exit program" << endl;
		return -1;
	}
	
	stringstream stream;
	getline(glb, text);
	getline(glb, text);
	stream << text;
	stream >> ndim >> nnode >> ngroup >>  nelem >> nmat>>nstep;
	stream.clear();
	glb.close();
	return 0;
}
int FEMOOP::ReadAllFile()
{
	ifstream cor(probn + ".cor");
	ifstream grp(probn + ".grp");
	ifstream ele(probn + ".ele");
	ifstream loa(probn + ".loa");
	ifstream pre(probn + ".pre");

	Nodes.ReadFile(cor,ndim,nnode);
	Groups.ReadFile(grp, ngroup);
	Elems.ReadFile(ele, Groups,ngroup); 
	Mats.ReadFile(probn,nmat);
	Loads = new Load[nstep];
	Pres = new Prescribe[nstep];
	for (int istep = 0; istep < nstep; istep++)
	{
		Loads[istep].ReadFile(loa);
		Pres[istep].ReadFile(pre);
	}
	return 0;
}
int FEMOOP::OpenCheckFile()
{
	chk.open(probn + ".chk");
	cout << "Project file located at " << probn << endl;
	chk << "Project file located at " << probn << endl;
	return 0;
}
int FEMOOP::ShowPassTime()
{
	time(&tNow);
	localtime_s(&tmLocal, &tNow);
	cout << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon + 1 << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0') << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
	chk << "Current Time " << tmLocal.tm_year + 1900 << "-" << tmLocal.tm_mon << "-" <<
		tmLocal.tm_mday << setw(4) << tmLocal.tm_hour << ":" << setfill('0')  << setw(2) << tmLocal.tm_min << ":"
		<< setfill('0') << setw(2) << tmLocal.tm_sec << setfill(' ') << endl;
	return 0;
}
int FEMOOP::TotalDOF(int istep,bool Rebuild)
{
	DegreeOfFreedom = new int[nnode*ndim]();
	int *Dof;
	Dof = new int[nnode*ndim];
	Rebuild = false;
	for (int idof = 0; idof < nnode*ndim; idof++)
	{
		Dof[idof] = 1;
	}
	TotalDegreeOfFreedom=Pres[istep].FixDof(Dof, nnode*ndim);
	for (int idof = 0; idof < nnode*ndim; idof++)
	{
		if (Dof[idof] != DegreeOfFreedom[idof])
		{
			DegreeOfFreedom[idof] = Dof[idof];
			Rebuild = true;
		}
	}
	delete[] Dof;
	Elems.FillDof(DegreeOfFreedom);
	return 0;
}
int FEMOOP::GlobalStiff()
{
	GlobalStiffMatrix = new double*[TotalDegreeOfFreedom];
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		GlobalStiffMatrix[iFreedom] = new double[TotalDegreeOfFreedom]();
	}
	int iedof,jedof,*Dof;
	double **Stiff;
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		Stiff = Elems.GetStiff(ielem);
		Dof = Elems.GetDof(ielem);
		for (int idof = 0; idof < Elems.Ndof(ielem); idof++)
		{
			iedof = Dof[idof];
			for (int jdof = 0; jdof < Elems.Ndof(ielem); jdof++)
			{
				jedof = Dof[jdof];
				if (iedof != 0 && jedof != 0)
				{
					GlobalStiffMatrix[iedof - 1][jedof - 1] += Stiff[idof][jdof];
				}
			}
		}
	}
	return 0;
}
int FEMOOP::ApplyLoad(int istep)
{
	//waiting to fill
	LoadMatrix[istep] = new double[TotalDegreeOfFreedom]();
	DispLoad = new double[TotalDegreeOfFreedom]();
	InteractLoad = new double[TotalDegreeOfFreedom]();
	InitialDisplace[istep] = new double[nnode*ndim]();
	Loads[istep].ApplyLoad(LoadMatrix[istep],DegreeOfFreedom);
	Pres[istep].ApplyDisp(DispLoad, DegreeOfFreedom, TotalDegreeOfFreedom, nnode, 2, nelem,InitialDisplace[istep]);

	return 0;
}
int FEMOOP::Solve(int istep)
{
	//waiting to fill
	int *ipiv, iteration;
	double *InitialForce0, *RightHand0,*InitialForce1;
	double *InteractDisplacement; 
	int nInteract, *InteractNode,*InteractNodeReceive ,InteractProcess,InteractnNode,tag1,tag2;
	double *InteractResult,*InteractReceive;
	char trana, tranb;
	int m, n, k;
	double alpha = 1.0, beta = 1.0;
	int lda, ldb, ldc;
	trana = 'N';
	tranb = 'N';
	InitialForce0 = new double[TotalDegreeOfFreedom]();
	InitialForce1 = new double[TotalDegreeOfFreedom]();
	ipiv = new int[TotalDegreeOfFreedom];
	RightHand0 = new double[TotalDegreeOfFreedom]();
	RightHand[istep] = new double[TotalDegreeOfFreedom]();
	InitialForce = new double[TotalDegreeOfFreedom]();
	InteractDisplacement = new double[TotalDegreeOfFreedom]();
	A = new double[TotalDegreeOfFreedom*TotalDegreeOfFreedom];
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		for (int jFreedom = 0; jFreedom < TotalDegreeOfFreedom; jFreedom++)
		{
			A[iFreedom*TotalDegreeOfFreedom + jFreedom] = GlobalStiffMatrix[iFreedom][jFreedom];
			//chk << setw(15) << GlobalStiffMatrix[iFreedom][jFreedom];
		}
		//chk << endl;
	}
	if (istep == 0)
	{
		Elems.ElementStress();
		Elems.InitialStrain(InitialForce);
		Elems.InitialStress(InitialForce);
	}
	iteration = 0;
	int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, TotalDegreeOfFreedom, TotalDegreeOfFreedom, A, TotalDegreeOfFreedom, ipiv);
	do
	{
		do
		{
			iteration++;
			cout << setw(20) << "istep" << setw(10) << istep << setw(20) << "iteration" << setw(10) << iteration << endl;
			chk << setw(20) << "istep" << setw(10) << istep << setw(20) << "iteration" << setw(10) << iteration << endl;
			for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
			{
				InitialForce0[iFreedom] = 0;
			}
			Nodes.SetZero("All");
			Elems.ElementLoad(InitialForce0);
			//cout << setw(15) << "LoadMatrix"
			//	<< setw(15) << "DispLoad"
			//	<< setw(15) << "InitialForce"
			//	<< setw(15) << "InitialForce0"
			//	<< setw(15) << "RightHand0" << endl;
			for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
			{
				RightHand0[iFreedom] = LoadMatrix[istep][iFreedom] + DispLoad[iFreedom] + InitialForce[iFreedom] - InitialForce0[iFreedom];
				//cout << setw(15) << LoadMatrix[istep][iFreedom]
				//	<< setw(15) << DispLoad[iFreedom]
				//	<< setw(15) << InitialForce[iFreedom]
				//	<< setw(15) << InitialForce0[iFreedom]
				//	<< setw(15) << RightHand0[iFreedom] << endl;
			}

			info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, trana, TotalDegreeOfFreedom, 1, A, TotalDegreeOfFreedom, ipiv, RightHand0, TotalDegreeOfFreedom);
			//int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, TotalDegreeOfFreedom, 1, A, TotalDegreeOfFreedom, ipiv, RightHand0, TotalDegreeOfFreedom);

			for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
			{
				RightHand[istep][iFreedom] += RightHand0[iFreedom];
			}
			Elems.GetResult(RightHand[istep]);
			Elems.GetInitialDisplacement(InitialDisplace[istep]);
			Nodes.SetZero("All");
			Elems.ElementStress();
			Elems.ElementLoad(InitialForce1);

			Error = 0;
			//cout << setw(15) << "LoadMatrix"
			//	<< setw(15) << "InitialForce"
			//	<< setw(15) << "DispLoad"
			//	<< setw(15) << "InitialForce1" << setw(10)<<myid<<endl;
			for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
			{
				//Error += pow((LoadMatrix[istep][iFreedom] + InitialForce[iFreedom] + DispLoad[iFreedom] + InitialForce0[iFreedom] - InitialForce1[iFreedom]), 2);
				Error += pow((LoadMatrix[istep][iFreedom] + InitialForce[iFreedom] + DispLoad[iFreedom] - InitialForce1[iFreedom]), 2);
				//cout << setw(15) << LoadMatrix[istep][iFreedom]
				//	<< setw(15) << InitialForce[iFreedom]
				//	<< setw(15) << DispLoad[iFreedom]
				//	<< setw(15) << InitialForce1[iFreedom] << endl;
				ErrorSum += pow((LoadMatrix[istep][iFreedom] + InitialForce[iFreedom] + DispLoad[iFreedom]), 2);
				InitialForce1[iFreedom] = 0;
			}
			ErrorSum = ErrorSum / TotalDegreeOfFreedom;
			Error = Error / ErrorSum;
			chk << setw(20) << "Load Error=" << setw(10) << Error << endl;
		} while (Error > 1e-11);
		if (istep > 0)
		{
			for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
			{
				RightHand[istep][iFreedom] += RightHand[istep - 1][iFreedom];
			}
		}
		nInteract = Pres[istep].GetnInteract();
		//for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
		//{
		//	DispLoad[iFreedom] -= InteractLoad[iFreedom];
		//}
		for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
		{
			InteractLoad[iFreedom] = 0;
		}
		cout << setw(20) << "Initial Displacement" <<setw(10)<< myid << endl;
		for (int iInteract = 0; iInteract < nInteract; iInteract++)
		{
			InteractnNode = Pres[istep].GetInteractnNode(iInteract);
			InteractProcess = Pres[istep].GetInteractProcess(iInteract);
			InteractNode = Pres[istep].GetInteractNode(iInteract);
			InteractNodeReceive = new int[InteractnNode]();
			InteractResult = new double[InteractnNode*ndim]();
			InteractReceive = new double[InteractnNode*ndim]();
			for (int inode = 0; inode < InteractnNode; inode++)
			{
				int NodeIndex = InteractNode[inode];
				for (int idim = 0; idim < ndim; idim++)
				{
				//int idim = 0;
					if (DegreeOfFreedom[NodeIndex * 2 + idim] != 0)
					{
						//cout << setw(20) << RightHand[istep][DegreeOfFreedom[NodeIndex * 2 + idim] - 1] ;
						InteractResult[inode * 2 + idim] = RightHand[istep][DegreeOfFreedom[NodeIndex * 2 + idim] - 1];
						/*RightHand[istep][DegreeOfFreedom[NodeIndex * 2 + idim] - 1] -= InteractResult[inode * 2 + idim];*/
						cout << setw(20) << RightHand[istep][DegreeOfFreedom[NodeIndex * 2 + idim] - 1] << setw(10) << myid << endl;
					}
					else
					{
						InteractResult[inode * 2 + idim] = 0;
					}
					//cout <<setw(15)<< InteractResult[inode * 2 +idim] << endl;
				}
			}
			tag1 = 1;
			tag2 = 2;
			//MPI::COMM_WORLD.Send(InteractResult, InteractnNode*ndim, MPI::DOUBLE, 0, tag1);
			//MPI::COMM_WORLD.Recv(InteractReceive, InteractnNode*ndim, MPI::DOUBLE, 0, tag1);
			MPI::COMM_WORLD.Send(InteractResult, InteractnNode*ndim, MPI::DOUBLE, InteractProcess, tag1); 
			MPI::COMM_WORLD.Recv(InteractReceive, InteractnNode*ndim, MPI::DOUBLE, InteractProcess, tag1);
			MPI::COMM_WORLD.Send(InteractNode, InteractnNode, MPI::INT, InteractProcess, tag2);
			MPI::COMM_WORLD.Recv(InteractNodeReceive, InteractnNode, MPI::INT, InteractProcess, tag2);
			//InteractNode = Pres[istep].GetNode(iInteract);
			//for (int inode = 0; inode < InteractnNode; inode++)
			//{
			//	if (InteractNode[inode] != InteractNodeReceive[inode])
			//	{
			//		cout << InteractNode[inode] << InteractNodeReceive[inode] << InteractReceive[inode*2] << endl;
			//		cout << -2 << endl;return -2;
			//	}
			//}
			for (int inode = 0; inode < InteractnNode; inode++)
			{
				//cout << InteractResult[inode] << setw(10)<<myid<<endl;
			}
			for (int inode = 0; inode < InteractnNode; inode++)
			{
				int NodeIndex = InteractNode[inode]; 
				for (int idim = 0; idim < ndim; idim++)
				{
					InteractResult[inode * 2 + idim] = (InteractResult[inode * 2 + idim] + InteractReceive[inode * 2 + idim]) / 2;
					//cout << InteractResult[inode * 2 + idim]<<endl;
				}
			}
			for (int inode = 0; inode < InteractnNode; inode++)
			{
				//cout << InteractResult[inode] << setw(10) << myid << endl;
			}


			cout << setw(20) << "Modify Displacement" << setw(10) << myid << endl;
			for (int inode = 0; inode < InteractnNode; inode++)
			{
				int NodeIndex = InteractNode[inode];
				for (int idim = 0; idim < ndim; idim++)
				{
					if (DegreeOfFreedom[NodeIndex * 2 + idim] != 0)
					{
						RightHand[istep][DegreeOfFreedom[NodeIndex * 2 + idim] - 1]=InteractResult[inode * 2 + idim];
						cout << setw(20) << RightHand[istep][DegreeOfFreedom[NodeIndex * 2 + idim] - 1] << setw(10) << myid << endl;
					}
				}
			}
			Pres[istep].ApplyInterDisp(InteractLoad, GlobalStiffMatrix, DegreeOfFreedom, TotalDegreeOfFreedom, nnode, 2, nelem, InteractDisplacement, InteractResult, InteractNode);
		}
		cout <<setw(10)<< "Dload" << endl;
		for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
		{
			DispLoad[iFreedom] += InteractLoad[iFreedom];
			cout << setw(10)<<InteractLoad[iFreedom]<<endl;
		}
		//for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
		//{
		//	  InitialDisplace[istep][iFreedom]+=InteractDisplacement[iFreedom];
		//}

		double *zero;
		zero=new double[TotalDegreeOfFreedom]();
		Elems.GetResult(RightHand[istep]);
		Elems.GetInitialDisplacement(InitialDisplace[istep]);
		Nodes.SetZero("All");
		Elems.ElementStress();
		Elems.ElementLoad(InitialForce1);
		Error = 0;
		cout << setw(15) << "LoadMatrix"
			<< setw(15) << "InitialForce"
			<< setw(15) << "DispLoad"
			<< setw(15) << "InitialForce1" << endl;
		for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
		{
			//Error += pow((LoadMatrix[istep][iFreedom] + InitialForce[iFreedom] + DispLoad[iFreedom] + InitialForce0[iFreedom] - InitialForce1[iFreedom]), 2);
			Error += pow((LoadMatrix[istep][iFreedom] + InitialForce[iFreedom] + DispLoad[iFreedom] - InitialForce1[iFreedom]), 2);
			cout << setw(15) << LoadMatrix[istep][iFreedom]
				<< setw(15) << InitialForce[iFreedom]
				<< setw(15) << DispLoad[iFreedom]
				<< setw(15) << InitialForce1[iFreedom] << endl;
			ErrorSum += pow((LoadMatrix[istep][iFreedom] + InitialForce[iFreedom] + DispLoad[iFreedom]), 2);
			InitialForce1[iFreedom] = 0;
		}
		ErrorSum = ErrorSum / TotalDegreeOfFreedom;
		Error = Error / ErrorSum;
		chk << setw(20) << "Load Error=" << setw(10) << Error << endl;
		cout << setw(20) << "Load Error=" << setw(10) << Error << endl;
		if (myid == 0)
		{
			char a;
			cin >> a;
			MPI::COMM_WORLD.Barrier();
		}
		else
		{
			MPI::COMM_WORLD.Barrier();
		}
	}while (Error > 1e-11);                    
	Elems.GetResult(RightHand[istep]);
	Elems.GetInitialDisplacement(InitialDisplace[istep]);
	Nodes.SetZero("All");
	Elems.ElementStress();
	Nodes.SetZero("Times");
	Elems.ElementLoad(InitialForce1);
	return 0;
}
int FEMOOP::PostProcess(int istep)
{
	//waiting to fill
	chk << setw(10) << "Index"
		<< setw(15) << "Dispx"
		<< setw(15) << "Dispy"
		<< setw(15) << "Dispz"
		<< endl; 
	double *disp;
	disp = new double[ndim];
	for (int inode = 0; inode < nnode; inode++)
	{		
		chk << setw(10) << inode;
		for (int idim = 0; idim < ndim; idim++)
		{
			if (DegreeOfFreedom[inode * 2 + idim]!=0)
			{
				disp[idim] = RightHand[istep][DegreeOfFreedom[inode * 2 + idim] - 1];
				chk << setw(15) << RightHand[istep][DegreeOfFreedom[inode * 2 + idim] - 1];
			}
			else
			{
				disp[idim] = 0;
				chk << setw(15) << 0;
			}
			disp[idim] += InitialDisplace[istep][inode * 2 + idim];
		}
		Nodes.PutResult(inode, disp,"Disp");
		chk << endl;
	}
	chk << setw(10) << "Index"
		<< setw(15) << "Dispx"
		<< setw(15) << "Dispy"
		<< setw(15) << "Dispz"
		<< endl;
	Nodes.PrintResult(chk);
	cout << setw(10) << "Index"
		<< setw(15) << "StressXX"
		<< setw(15) << "StressYY"
		<< setw(15) << "StressXY" << endl;
	chk << setw(10) << "Index"
		<< setw(15) << "StressXX"
		<< setw(15) << "StressYY"
		<< setw(15) << "StressXY" << endl;
	
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		Elems.PrintStress(ielem, chk);
		Elems.PrintStress(ielem, cout);
	}
	return 0;
}
int FEMOOP::StrainForce()
{
	return 0;
}
int FEMOOP::GIDOutMesh()
{
	string GidMeshFile;
	double *Coor;
	int *ENode;
	GidMeshFile = probn + ".flavia.msh";
	GiD_OpenPostMeshFile(GidMeshFile.c_str(), GiD_PostAscii);
	GiD_BeginMesh("Static", GiD_2D, GiD_Quadrilateral, 4);
	GiD_BeginCoordinates();
	
	for (int inode = 0; inode < nnode; inode++)
	{
		Coor=Nodes.GetCoor(inode);
		GiD_WriteCoordinates2D(inode + 1, Coor[0], Coor[1]);
	}
	GiD_EndCoordinates();
	GiD_BeginElements();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		int Type = Elems.GetType(ielem);
		if (Type == 4)
		{
			ENode = new int[4];
			for (int inode = 0; inode < 4; inode++)
			{
				ENode[inode] = Elems.GetNode(ielem, inode)+1;
			}
			GiD_WriteElement(ielem + 1, ENode);
		}
		
	}
	GiD_EndElements();
	GiD_EndMesh();
	GiD_ClosePostMeshFile();

	return 0;
}
int FEMOOP::GIDOutResult(int istep)
{
	double *Stress, *Strain, *Disp,*InternalForce;
	
	GiD_BeginResult("Displacement", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nnode; inode++)
	{
		Disp = Nodes.GetDisp(inode);
		GiD_Write2DVector(inode + 1, Disp[0], Disp[1]);
	}
	GiD_EndResult();
	GiD_BeginResult("InternalForce", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nnode; inode++)
	{
		InternalForce = Nodes.GetInternalForce(inode);
		GiD_Write2DVector(inode + 1, InternalForce[0], InternalForce[1]);
	}
	GiD_EndResult();
	GiD_BeginResult("Strain", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nnode; inode++)
	{
		Strain = Nodes.GetStrain(inode);
		GiD_Write2DMatrix(inode + 1, Strain[0], Strain[1], Strain[2]);
	}
	GiD_EndResult();
	GiD_BeginResult("Stress", "Static", istep, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL);
	for (int inode = 0; inode < nnode; inode++)
	{
		Stress = Nodes.GetStress(inode);
		GiD_Write2DMatrix(inode + 1, Stress[0], Stress[1], Stress[2]);
	}
	GiD_EndResult();

	return 0;
}
int FEMOOP::StepCycle()
{
	//InputFile();
	if (!MPI::Is_initialized())
	{
		return -1; 
	}

	OpenCheckFile();
	ShowPassTime();
	ReadControl();
	ReadAllFile();
	ShowPassTime();
	GIDOutMesh();
	
	LoadMatrix = new double*[nstep];
	RightHand = new double *[nstep];
	InitialDisplace = new double*[nstep];
	GidResFile = probn + ".flavia.res";
	GiD_OpenPostResultFile(GidResFile.c_str(), GiD_PostAscii);
	TotalDOF(0, Rebuild);
	Elems.ElementStiff();
	GlobalStiff();
	ShowPassTime();
	MPI::COMM_WORLD.Barrier();

	for (int istep = 0; istep < nstep; istep++)
	{
		
		TotalDOF(istep,Rebuild);
		if (Rebuild)
		{
			delete[] GlobalStiffMatrix;
			Elems.ElementStiff();
			GlobalStiff();
			ShowPassTime();
		}		

		ApplyLoad(istep);
		ShowPassTime();

		Solve(istep);
		ShowPassTime();

		PostProcess(istep);
		
		GIDOutResult(istep);
		
		ShowPassTime();
	}
	GiD_ClosePostResultFile();
	ShowPassTime();
	return 0;
}