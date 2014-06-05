#include "Prescribe.h"
#include "mkl.h"
#include <algorithm>
Prescribe *Pres;
extern Element Elems;
Prescribe::Prescribe()
{
}
FixBoundary::FixBoundary(string text)
{
	stringstream stream;
	stream << text;
	stream >> Dir >> Nnode;
	Dir--;
	stream.clear();
	Node = new int[Nnode];
}
int FixBoundary::AddVal(string text)
{
	stringstream stream;
	stream << text;
	for (int inode = 0; inode < Nnode; inode++)
	{
		stream >> Node[inode];
		Node[inode]--;
	}
	return 0;
}
Displacement::Displacement(string text)
{
	stringstream stream;
	stream << text;
	stream >> Dir >> Value >> Nnode;
	Dir--;
	stream.clear();
	Node = new int[Nnode];
}
int Displacement::AddNode(string text)
{
	stringstream stream;
	stream << text;
	for (int inode = 0; inode < Nnode; inode++)
	{
		stream >> Node[inode];
		Node[inode]--;
	}
	return 0;
}
InteractBoundary::InteractBoundary(string text)
{
	stringstream stream;
	stream << text;
	stream >> nNode>>InteractDomain;
	Node = new int[nNode];
	InteractNode = new int[nNode];
	InteractDomain--;
}
int InteractBoundary::AddNode(string text)
{
	stringstream stream;
	stream << text;
	for (int inode = 0; inode < nNode; inode++)
	{
		stream >> Node[inode];
		Node[inode]--; 
	}
	return 0;
}
int InteractBoundary::AddInteractNode(string text)
{
	stringstream stream;
	stream << text;
	for (int inode = 0; inode < nNode; inode++)
	{
		stream >> InteractNode[inode];
		InteractNode[inode]--;
	}
	return 0;
}
int Prescribe::ReadFile(ifstream &pre)
{
	
	string text;
	stringstream stream;
	int iFix, iDisp,iInteract;
	getline(pre, text);
	getline(pre, text);
	stream << text;
	stream >> npre;
	stream.clear();
	iFix = 0; iDisp = 0, iInteract=0;
	for (int ipre = 0; ipre < npre; ipre++)
	{
		getline(pre, text);
		type.push_back(text);
		if (text == "FixBoundary")
		{
			getline(pre, text);
			stream << text;
			stream >> ndim;
			stream.clear();
			Index.push_back(iFix);
			iFix++;
			for (int idim = 0; idim < ndim; idim++)
			{
				getline(pre, text);
				getline(pre, text);
				stream << text;
				Fix.push_back(text);
				getline(pre, text);
				Fix[idim].AddVal(text);
			}
			stream.clear();
		}
		else if (text == "Displacement")
		{
			getline(pre, text);
			stream << text;
			stream >> nDisp;
			stream.clear();
			Index.push_back(iDisp);
			iDisp++;
			for (int idisp = 0; idisp < nDisp; idisp++)
			{
				getline(pre, text);
				getline(pre, text);
				stream << text;
				Disp.push_back(text);
				getline(pre, text);
				Disp[idisp].AddNode(text);
				stream.clear();
			}
		}
		else if (text == "InteractBoundary")
		{
			stream.str("");
			getline(pre, text);
			stream << text;
			stream >> nInteract;
			Index.push_back(iInteract);
			iInteract++;
			for (int iinteract = 0; iinteract < nInteract; iinteract++)
			{
				getline(pre, text);
				getline(pre, text);
				stream << text;
				Interact.push_back(text);
				getline(pre, text);
				Interact[iinteract].AddNode(text);
				stream.clear();
				getline(pre, text);
				Interact[iinteract].AddInteractNode(text);
			}
		}
	}
	return 0;
}
int FixBoundary::Get(string text)
{
	if (text == "Dir")
	{
		return Dir;
	}
	else if (text == "Nnode")
	{
		return Nnode;
	}
	return -1;
}
int FixBoundary::Get(int index)
{
	return Node[index];
}
int Displacement::GetNnode()
{
	return Nnode;
}
int Displacement::GetNode(int inode)
{
	return Node[inode];
}
int Displacement::GetDir()
{
	return Dir;
}
int InteractBoundary::GetNnode()
{
	return nNode;
}
int InteractBoundary::GetNode(int inode)
{
	return Node[inode];
}
int InteractBoundary::GetInteractNode(int inode)
{
	return InteractNode[inode];
}
int Prescribe::FixDof(int *DegreeOfFreedom, int n)
{
	for (int ipre = 0; ipre < npre; ipre++)
	{
		if (type.at(ipre) == "FixBoundary")
		{
			for (int idim = 0; idim < ndim; idim++)
			{
				if (idim == Fix[idim].Get("Dir") )
				{
					for (int inode = 0; inode < Fix[idim].Get("Nnode"); inode++)
					{
						DegreeOfFreedom[(Fix[idim].Get(inode)) * 2 + idim] = 0;
					}
				}
			}
		}
		//else if (type.at(ipre) == "InteractBoundary")
		//{
		//	int nBoundary = Interact.size();
		//	for (int iBoundary = 0; iBoundary < nBoundary; iBoundary++)
		//	{
		//		int nnode = Interact[iBoundary].GetNnode();
		//		for (int inode = 0; inode < nnode; inode++)
		//		{
		//			int node = Interact[iBoundary].GetNode(inode);
		//			DegreeOfFreedom[node * 2] = 0;
		//			DegreeOfFreedom[node * 2 + 1] = 0;
		//		}
		//	}
		//}
		else if (type.at(ipre) == "Displacement")
		{
			int nnode = Disp[Index[ipre]].GetNnode();
			int NodeIndex,iFreedom,Dir;
			Dir = Disp[Index[ipre]].GetDir();
			for (int inode = 0; inode < nnode; inode++)
			{
				NodeIndex = Disp[Index[ipre]].GetNode(inode);
				iFreedom = NodeIndex * 2 + Dir;
				DegreeOfFreedom[iFreedom] = 0;
			}
		}

	}
	int TotalDegreeOfFreedom = 0;
	for (int idof = 0; idof < n; idof++)
	{
		if (DegreeOfFreedom[idof] != 0)
		{
			TotalDegreeOfFreedom++;
			DegreeOfFreedom[idof] = TotalDegreeOfFreedom;
		}
	}

	return TotalDegreeOfFreedom;
}
int Prescribe::ApplyDisp(double *LoadMatrix, int *DegreeOfFreedom, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement)
{
	for (int idisp = 0; idisp < nDisp; idisp++)
	{
		Disp[idisp].ApplyDisp(LoadMatrix, DegreeOfFreedom, TotalDegreeOfFreedom, TotalNode, NodeDof,nelem,InitialDisplacement);
	}
	return 0;
}
int Displacement::ApplyDisp(double *LoadMatrix, int *DegreeOfFreedom, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement)
{
	double *NodeDisplacement,*DispLoad, *A;
	double **StiffMatrix;
	int NodeIndex,iDof,nDispDof,*DispDof;
	NodeDisplacement = new double[Nnode]();
	DispLoad = new double[TotalDegreeOfFreedom]();
	A = new double[TotalDegreeOfFreedom*TotalDegreeOfFreedom];

	int nDofMatrix = TotalNode*NodeDof;
	DispDof = new int[nDofMatrix]();
	for (int iFreedom = 0; iFreedom < nDofMatrix; iFreedom++)
	{
		DispDof[iFreedom] = DegreeOfFreedom[iFreedom];
	}
	nDispDof = TotalDegreeOfFreedom;
	for (int inode = 0; inode < Nnode; inode++)
	{
		NodeIndex = Node[inode];
		nDispDof++;
		DispDof[NodeIndex * 2 + Dir] = nDispDof;
		NodeDisplacement[inode] = Value;
		InitialDisplacement[NodeIndex * 2 + Dir] = Value;
	}
	Elems.FillDof(DispDof);
	nDispDof -= TotalDegreeOfFreedom;
	int iedof, jedof, *Dof;
	double **Stiff;
	StiffMatrix = new double*[TotalDegreeOfFreedom];
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		StiffMatrix[iFreedom] = new double[nDispDof]();
	}
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
					if (jedof>TotalDegreeOfFreedom &&iedof<=TotalDegreeOfFreedom)
					{
						StiffMatrix[iedof - 1][jedof - 1-TotalDegreeOfFreedom] += Stiff[idof][jdof];
					}
				}
			}
		}
	}
	Elems.FillDof(DegreeOfFreedom);
	double a, b, c;
	char trana, tranb;
	int m, n, k;
	double alpha = 1.0, beta = 1.0;
	trana = 'N';
	tranb = 'N';
	int lda, ldb, ldc;
	m = TotalDegreeOfFreedom; n = 1; k = nDispDof;
	lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		for (int jFreedom = 0; jFreedom < nDispDof; jFreedom++)
		{
			A[jFreedom*TotalDegreeOfFreedom + iFreedom] = StiffMatrix[iFreedom][jFreedom];
			//cout << setw(15) << StiffMatrix[iFreedom][jFreedom];
		}
		//cout << endl;
	}

	dgemm(&trana, &tranb, &m, &n, &k, &alpha, A, &lda, NodeDisplacement, &ldb, &beta, DispLoad, &ldc);
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		LoadMatrix[iFreedom] -= DispLoad[iFreedom];
	}
	
	return 0;
}
int Prescribe::ApplyInterDisp(double *LoadMatrix, double  **StiffMatrix, int *DegreeOfFreedom, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement, double *InteractResult, int *Interactnode)
{
	for (int iInteract = 0; iInteract < nInteract; iInteract++)
	{
		Interact[iInteract].ApplyInterDisp(LoadMatrix, StiffMatrix,DegreeOfFreedom, TotalDegreeOfFreedom, TotalNode, NodeDof, nelem, InitialDisplacement, InteractResult, Interactnode);
	}
	return 0;
}
int InteractBoundary::ApplyInterDisp(double *LoadMatrix, double  **StiffMatrix,int *DegreeOfFreedom, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement, double *InteractResult, int *Interactnode)
{
	double *NodeDisplacement, *DispLoad, *A;
	double **DispStiffMatrix;
	int NodeIndex, iDof, nDispDof, *DispDof;
	double *ForceDisp;
	NodeDisplacement = new double[nNode*2]();
	DispLoad = new double[TotalDegreeOfFreedom]();
	A = new double[TotalDegreeOfFreedom*TotalDegreeOfFreedom];
	DispStiffMatrix = new double*[TotalDegreeOfFreedom];
	ForceDisp = new double[TotalDegreeOfFreedom]();
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		DispStiffMatrix[iFreedom] = new double[TotalDegreeOfFreedom]();
	}
	int nDofMatrix = TotalNode*NodeDof;
	//DispDof = new int[nDofMatrix]();
	//for (int iFreedom = 0; iFreedom < nDofMatrix; iFreedom++)
	//{
	//	DispDof[iFreedom] = DegreeOfFreedom[iFreedom];
	//}
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		InitialDisplacement[iFreedom] = 0;
	}
	for (int inode = 0; inode < nNode; inode++)
	{
		NodeIndex = Node[inode];
		for (int idim = 0; idim < 2; idim++)
		{
			
			//DispDof[NodeIndex * 2 + idim] = nDispDof;
			NodeDisplacement[inode * 2 + idim] = InteractResult[inode * 2 + idim];
			iDof = DegreeOfFreedom[NodeIndex * 2 + idim]-1;
			InitialDisplacement[iDof] = InteractResult[inode * 2 + idim];
		}
	}
	//Elems.FillDof(DispDof);
	//nDispDof -= TotalDegreeOfFreedom;
	int iedof, jedof, *Dof;
	/*double **Stiff;
	StiffMatrix = new double*[TotalDegreeOfFreedom];
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		StiffMatrix[iFreedom] = new double[nDispDof]();
	}
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
					if (jedof>TotalDegreeOfFreedom &&iedof <= TotalDegreeOfFreedom)
					{
						StiffMatrix[iedof - 1][jedof - 1 - TotalDegreeOfFreedom] += Stiff[idof][jdof];
					}
				}
			}
		}
	}
	Elems.FillDof(DegreeOfFreedom);*/
	double a, b, c;
	char trana, tranb;
	int m, n, k;
	int *ipiv;
	double alpha = 1.0, beta = 1.0;
	trana = 'N';
	tranb = 'N';
	int lda, ldb, ldc;
	m = TotalDegreeOfFreedom; n = 1; k = TotalDegreeOfFreedom;
	lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
	ipiv = new int[TotalDegreeOfFreedom];
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		for (int jFreedom = 0; jFreedom < TotalDegreeOfFreedom; jFreedom++)
		{
			A[jFreedom*TotalDegreeOfFreedom + iFreedom] = StiffMatrix[iFreedom][jFreedom];
			//cout << setw(12) << StiffMatrix[iFreedom][jFreedom];
		}
		//cout << endl;
	}

	dgemm(&trana, &tranb, &m, &n, &k, &alpha, A, &lda, InitialDisplacement, &ldb, &beta, DispLoad, &ldc);
	nDispDof = 0;
	DispDof = new int[TotalDegreeOfFreedom]();
	for (int inode = 0; inode < nNode; inode++)
	{
		NodeIndex = Node[inode];
		for (int idim = 0; idim < 2; idim++)
		{
			iDof = DegreeOfFreedom[NodeIndex * 2 + idim] - 1;
			if (InitialDisplacement[iDof] != 0)
			{
				nDispDof++;
				DispDof[iDof] = 1;
				DispLoad[iDof] = -DispLoad[iDof];
			}
		}
	}
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		DispLoad[iFreedom] = -DispLoad[iFreedom];
	}
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		if (DispDof[iFreedom] != 0)
		{
			for (int jFreedom = 0; jFreedom < TotalDegreeOfFreedom; jFreedom++)
			{
				DispStiffMatrix[iFreedom][jFreedom] = 0;
				if (iFreedom == jFreedom)
				{
					DispStiffMatrix[iFreedom][jFreedom] = 1;
				}
			}
		}
		else
		{
			for (int jFreedom = 0; jFreedom < TotalDegreeOfFreedom; jFreedom++)
			{
				DispStiffMatrix[iFreedom][jFreedom] = StiffMatrix[iFreedom][jFreedom];
				if (DispDof[jFreedom] != 0)
				{
					DispStiffMatrix[iFreedom][jFreedom] = 0;
				}
			}
		}
	}
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		for (int jFreedom = 0; jFreedom < TotalDegreeOfFreedom; jFreedom++)
		{
			A[jFreedom*TotalDegreeOfFreedom + iFreedom] = DispStiffMatrix[iFreedom][jFreedom];
			//cout << setw(12) << StiffMatrix[iFreedom][jFreedom];
		}
		//cout << endl;
	}
	for (int iFreedom = 0; iFreedom< TotalDegreeOfFreedom; iFreedom++)
	{
		ForceDisp[iFreedom] = DispLoad[iFreedom];
	}
	int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, TotalDegreeOfFreedom, 1, A, TotalDegreeOfFreedom, ipiv, ForceDisp, TotalDegreeOfFreedom);
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		if (InitialDisplacement[iFreedom] == 0)
		{
			InitialDisplacement[iFreedom] = ForceDisp[iFreedom];
		}
	}
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		DispLoad[iFreedom] = 0;
	}
	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		for (int jFreedom = 0; jFreedom < TotalDegreeOfFreedom; jFreedom++)
		{
			A[jFreedom*TotalDegreeOfFreedom + iFreedom] = StiffMatrix[iFreedom][jFreedom];
			//cout << setw(12) << StiffMatrix[iFreedom][jFreedom];
		}
		//cout << endl;
	}

	dgemm(&trana, &tranb, &m, &n, &k, &alpha, A, &lda, InitialDisplacement, &ldb, &beta, DispLoad, &ldc);


	for (int iFreedom = 0; iFreedom < TotalDegreeOfFreedom; iFreedom++)
	{
		LoadMatrix[iFreedom] -= DispLoad[iFreedom];
	}

	return 0;
}
int Prescribe::GetnInteract()
{
	return nInteract;
}
int InteractBoundary::GetInteractDomain()
{
	return InteractDomain;
}
int Prescribe::GetInteractProcess(int iInteract)
{
	return Interact[iInteract].GetInteractDomain();
}
int Prescribe::GetInteractnNode(int iInteract)
{
	return Interact[iInteract].GetNnode();
}
int *Prescribe::GetInteractNode(int iInteract)
{
	int *Node,nNode;
	nNode = Interact[iInteract].GetNnode();

	Node = new int[nNode];
	for (int inode = 0; inode < nNode; inode++)
	{
		Node[inode] = Interact[iInteract].GetNode(inode);
	}
	return Node;
}
Prescribe::~Prescribe()
{
}
