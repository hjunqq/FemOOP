#include "Element.h"
#include "mkl.h"
#include <algorithm>
extern Node Nodes;
extern Group Groups;
extern Material Mats;
Element Elems;
const double Rstg[] = { -0.774596669241483, 0.0, 0.774596669241483 };
const double H[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };

Element::Element()
{

}
int Quadrilateral::EGetCoor()
{
	Coor = new double *[Nnode];
	for (int inode = 0; inode < Nnode; inode++)
	{
		Coor[inode] = Nodes.GetCoor(Node[inode]);
	}
	return 0;
}
int Quadrilateral::Shape(double xi, double eta)
{
	PartialDeri[0][0] = -0.25*(1 - eta);
	PartialDeri[1][0] = -0.25*(1 - xi);
	PartialDeri[0][1] = 0.25*(1 - eta);
	PartialDeri[1][1] = -0.25*(1 + xi);
	PartialDeri[0][2] = 0.25*(1 + eta);
	PartialDeri[1][2] = 0.25*(1 + xi);
	PartialDeri[0][3] = -0.25*(1 + eta);
	PartialDeri[1][3] = 0.25*(1 - xi);
	Function[0] = 0.25*(1 - xi)*(1 - eta);
	Function[1] = 0.25*(1 + xi)*(1 - eta);
	Function[2] = 0.25*(1 + xi)*(1 + eta);
	Function[3] = 0.25*(1 - xi)*(1 + eta);
	Det = 0.0;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Xjac[i][j] = 0.0;
			for (int k = 0; k < 4; k++)
			{
				Xjac[i][j] += PartialDeri[i][k] * Coor[k][j];
			}
		}
	}
	Det = Xjac[0][0] * Xjac[1][1] - Xjac[1][0] * Xjac[0][1];
	Rjac[0][0] = Xjac[1][1] / Det;
	Rjac[1][1] = Xjac[0][0] / Det;
	Rjac[1][0] = -Xjac[1][0] / Det;
	Rjac[0][1] = -Xjac[0][1] / Det;
	for (int k = 0; k < 4; k++)
	{
		for (int i = 0; i < 2; i++)
		{
			PDeri[i][k] = 0.0;
			for (int j = 0; j < 2; j++)
			{
				PDeri[i][k] = PDeri[i][k] + Rjac[i][j] * PartialDeri[j][k];
			}
		}
	}
	return 0;
}
int Quadrilateral::ElementStiff()
{
	double Elastic, Posion, Density;
	double DNIX, DNIY, DNJX, DNJY, Dxx, Dxy, Dyx, Dyy;
	Material = Groups.Get(Group, "Material");
	Elastic = Mats.Get(Material, "Modulus");
	Posion = Mats.Get(Material, "Posion");
	Density = Mats.Get(Material, "Density");
	double D1 = Elastic*(1.000 - Posion) / ((1.00 + Posion)*(1 - 2.00*Posion));
	double D2 = Elastic*Posion / ((1.00 + Posion)*(1.00 - 2.00*Posion));
	double D3 = Elastic*0.50 / (1.00 + Posion);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Dxx = Dxy = Dyx = Dyy = 0;
			for (int ieta = 0; ieta < 3; ieta++)
			{
				double eta = Rstg[ieta];
				double etah = H[ieta];
				for (int ixi = 0; ixi < 3; ixi++)
				{
					double xi = Rstg[ixi];
					double xih = H[ixi];
					Shape(xi, eta);
					DNIX = PDeri[0][i];
					DNIY = PDeri[1][i];
					DNJX = PDeri[0][j];
					DNJY = PDeri[1][j];
					Dxx = Dxx + DNIX * DNJX * Det * etah * xih;
					Dxy = Dxy + DNIX * DNJY * Det * etah * xih;
					Dyx = Dyx + DNIY * DNJX * Det * etah * xih;
					Dyy = Dyy + DNIY * DNJY * Det * etah * xih;
				}
			}
			Stiff[i * 2][j * 2] = Dxx*D1 + Dyy*D3;
			Stiff[i * 2 + 1][j * 2 + 1] = Dyy*D1 + Dxx*D3;
			Stiff[i * 2][j * 2 + 1] = Dxy*D2 + Dyx*D3;
			Stiff[i * 2 + 1][j * 2] = Dyx*D2 + Dxy*D3;
		}
	}
	return 0;
}
Quadrilateral::Quadrilateral(string text)
{
	stringstream stream;
	stream << text;
	stream >> Index >>
		Node[0] >>
		Node[1] >>
		Node[2] >>
		Node[3] >> Group;
	Index--;
	Node[0]--;
	Node[1]--;
	Node[2]--;
	Node[3]--;
	Group--;
	for (int i = 0; i < 3; i++)
	{
		Strain[i] = 0;
		Stress[i] = 0;
	}
	for (int i = 0; i < 8; i++)
	{
		Result[i] = 0;
	}
	Dof = new int[8]();
}
int Element::ReadFile(ifstream &ele, Group &Groups,int ngroup)
{
	string text;
	stringstream stream;
	int iQuadrs=0;
	for (int igroup = 0; igroup < ngroup; igroup++)
	{
		if (Groups.Get(igroup,"Type") == 4)
		{			
			for (int ielem = 0; ielem < Groups.Get(igroup,"Nelem"); ielem++)
			{
				Type.push_back(4);
				Groups.PutElem(igroup,ielem, Type.size() - 1);
				Index.push_back(iQuadrs);
				getline(ele, text);
				Quadrs.push_back(text);
				Quadrs[iQuadrs].EGetCoor();
				iQuadrs++;
			}	
		}
	}
	return 0;
}
int Element::ElementStiff()
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type.at(ielem) == 4)
		{
			Quadrs[Index[ielem]].ElementStiff();
		}
	}
	return 0;
}
int Element::ElementStress()
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type.at(ielem) == 4)
		{
			Quadrs[Index[ielem]].ElementStress();
		}
	}
	return 0;
}
double ** Quadrilateral::GetStiff()
{
	double **Stiff;
	Stiff = new double *[8];
	for (int i = 0; i < 8; i++)
	{
		Stiff[i] = new double[8];
		for (int j = 0; j < 8; j++)
		{
			Stiff[i][j] = this->Stiff[i][j];
		}
	}
	return Stiff;
}
int  * Quadrilateral::GetDof()
{

	return Dof;
}
double ** Element::GetStiff(int ielem)
{
	if (Type[ielem] == 4)
	{
		return Quadrs[Index[ielem]].GetStiff();
	}
	return 0;
}
int *Element::GetDof(int ielem)
{
	if (Type[ielem] == 4)
	{
		return Quadrs[Index[ielem]].GetDof();
	}
	return 0;
}
int Quadrilateral::GetNode(int inode)
{
	return Node[inode];
}
int Element::GetNode(int ielem, int inode)
{
	if (Type[ielem] == 4)
	{
		return Quadrs[Index[ielem]].GetNode(inode);
	}
	return 0;
}
int Quadrilateral::FillDof(int *DegreeOfFreedom)
{
	for (int inode = 0; inode < Nnode; inode++)
	{
		int index = Node[inode];
		for (int idim = 0; idim < 2; idim++)
		{
			Dof[inode * 2 + idim] = DegreeOfFreedom[(index ) * 2 + idim];
		}
	}
	return 0;
}
int Element::FillDof(int *DegreeOfFreedom)
{

	int	nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].FillDof(DegreeOfFreedom);
		}
	}
	return 0;
}
int Element::Ndof(int ielem)
{
	if (Type[ielem] == 4)
	{
		return 8;
	}
	return 0;
}
int Quadrilateral::GetResult(double *Result)
{
	for (int i = 0; i < 8; i++)
	{
		//this->Result[i] = 0;
		if (Dof[i] != 0)
		{
			this->Result[i] = Result[Dof[i]-1];
		}
	}
	return 0;
}
int Element::GetResult(double *Result)
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].GetResult(Result);
		}
	}
	return 0;
}
int Quadrilateral::GetInitialDisplacement(double *InitialDisplacement)
{
	for (int inode = 0; inode < 4; inode++)
	{
		for (int idim = 0; idim < 2; idim++)
		{
			Result[inode * 2 + idim] += InitialDisplacement[Node[inode] * 2 + idim];
		}
	}
	return 0;
}
int Quadrilateral::GetInteractDisp(double *InteractDisp, int *InteractNode, int nInteractNode)
{
	for (int inode = 0; inode < 4; inode++)
	{
		for (int idim = 0; idim < 2; idim++)
		{
			for (int iInterNode = 0; iInterNode < nInteractNode; iInterNode++)
			{
				if (Node[inode] == InteractNode[iInterNode])
				{
					for (int idim = 0; idim < 2; idim++)
					{
						Result[inode * 2+idim] = InteractDisp[iInterNode * 2+idim];
					}
				}
			}
		}
	}
	return 0;
}
int Element::GetInitialDisplacement(double *InitialDisplacement)
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].GetInitialDisplacement(InitialDisplacement);
		}
	}
	return 0;
}
int Element::GetInteractDisp(double *InteractDisp, int *InteractNode, int nInteractNode)
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].GetInteractDisp(InteractDisp,InteractNode,nInteractNode);
		}
	}
	return 0;
}
int Quadrilateral::ElementStress()
{
	double Elastic, Posion, Density;
	Material = Groups.Get(Group, "Material");
	Elastic = Mats.Get(Material, "Modulus");
	Posion = Mats.Get(Material, "Posion");
	Density = Mats.Get(Material, "Density");
	double D1 = Elastic*(1.000 - Posion) / ((1.00 + Posion)*(1 - 2.00*Posion));
	double D2 = Elastic*Posion / ((1.00 + Posion)*(1.00 - 2.00*Posion));
	double D3 = Elastic*0.50 / (1.00 + Posion);
	double eta = Rstg[1];
	double xi = Rstg[1];
	double xicor[] = { -0.577350269, 0.577350269, 0.577350269, -0.577350269 };
	double etacor[] = { -0.577350269, -0.577350269, 0.577350269, 0.577350269 };
	Shape(xi, eta);
	for (int i = 0; i < 3; i++)
	{
		Strain[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			NodeStrain[j][i] = 0;
			NodeStress[j][i] = 0;
			GaussStrain[j][i] = 0;
			GaussStress[j][i] = 0;
		}
	}
	for (int i = 0; i < 4; i++)
	{
		int ii = 2 * i;
		int j1 = ii;
		int j2 = ii + 1;
		double Bi = PDeri[0][i];
		double Ci = PDeri[1][i];
		B[0][j1] = Bi;
		B[1][j1] = 0;
		B[2][j1] = Ci;
		B[0][j2] = 0;
		B[1][j2] = Ci;
		B[2][j2] = Bi;
	}

	for (int k = 0; k < 4; k++)
	{
		for (int l = 0; l < 3; l++)
		{
			Strain[l] += B[l][2 * k] * Result[2 * k] + B[l][2 * k + 1] * Result[2 * k + 1];
		}
	}
	Stress[0] = D1*Strain[0] + D2*Strain[1];
	Stress[1] = D2*Strain[0] + D1*Strain[1];
	Stress[2] = D3*Strain[2];
	for (int i = 0; i < 4; i++)
	{
		xi = xicor[i];
		eta = etacor[i];
		Shape(xi, eta);
		for (int j = 0; j < 4; j++)
		{
			int ii = 2 * j;
			int j1 = ii;
			int j2 = ii + 1;
			double Bi = PDeri[0][j];
			double Ci = PDeri[1][j];
			B[0][j1] = Bi;
			B[1][j1] = 0;
			B[2][j1] = Ci;
			B[0][j2] = 0;
			B[1][j2] = Ci;
			B[2][j2] = Bi;
			for (int k = 0; k < 3; k++)
			{
				GaussStrain[i][k] += B[k][2 * j] * Result[2 * j] + B[k][2 * j + 1] * Result[2 * j + 1];
			}
		}
		GaussStress[i][0] = D1*GaussStrain[i][0] + D2*GaussStrain[i][1];
		GaussStress[i][1] = D2*GaussStrain[i][0] + D1*GaussStrain[i][1];
		GaussStress[i][2] = D3*GaussStrain[i][2];
	}

	//应力外推
	double a, b, c, *NValue0,*NValue1;
	char trana, tranb;
	int m, n, k;
	double alpha = 1.0, beta = 1.0;
	trana = 'N';
	tranb = 'N';
	int lda, ldb, ldc;
	a = 1 + sqrt(3) / 2; b = -1.0 / 2; c = 1 - sqrt(3) / 2;
	double TransMatrix[4][4] = { a, b, c, b,
		b, a, b, c,
		c, b, a, b,
		b, c, b, a };
	m = 4; n = 1; k = 4;
	lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
	NValue0 = new double[4];
	NValue1 = new double[4]();
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			NValue0[j] = GaussStress[j][i];
			NValue1[j] = 0;
		}
		dgemm(&trana, &tranb, &m, &n, &k, &alpha, TransMatrix[0], &lda, NValue0, &ldb, &beta, NValue1, &ldc);
		for (int j = 0; j < 4; j++)
		{
			NodeStress[j][i] = NValue1[j];
		}


		for (int j = 0; j < 4; j++)
		{
			NValue0[j] = GaussStrain[j][i];
			NValue1[j] = 0;
		}
		dgemm(&trana, &tranb, &m, &n, &k, &alpha, TransMatrix[0], &lda, NValue0, &ldb, &beta, NValue1, &ldc);
		for (int j = 0; j < 4; j++)
		{
			NodeStrain[j][i] = NValue1[j];
		}
		
	}
	double *disp;
	disp = new double[2];
	for (int i = 0; i < 4; i++)
	{
		disp[0] = Result[i*2];
		disp[1] = Result[i * 2 + 1];
		Nodes.PutResult(Node[i], NULL, "Time");
		Nodes.PutResult(Node[i], NodeStrain[i], "Strain");
		Nodes.PutResult(Node[i], NodeStress[i], "Stress");
		Nodes.PutResult(Node[i], disp, "Disp");
	}
	return 0;
}
int Quadrilateral::PrintStress(ofstream & outdist)
{
	outdist << setw(10) << Index
		<<setw(15)<< Stress[0] 
		<< setw(15) << Stress[1] 
		<< setw(15) << Stress[2] << endl;
	return 0;
}
int Quadrilateral::PrintStress(ostream & outdist)
{
	outdist << setw(10) << Index
		<< setw(15) << Stress[0]
		<< setw(15) << Stress[1]
		<< setw(15) << Stress[2] << endl;
	return 0;
}
int Element::PrintStress(int ielem, ofstream & outdist)
{

	if (Type[ielem] == 4)
	{
		Quadrs[Index[ielem]].PrintStress(outdist);
	}

	return 0;
}
int Element::PrintStress(int ielem, ostream & outdist)
{

	if (Type[ielem] == 4)
	{
		Quadrs[Index[ielem]].PrintStress(outdist);
	}

	return 0;
}
int Quadrilateral::InitialStress(double *InitialForce)
{
	double B[3][2], S[2][3] = { 0 };
	double *NStress;
	char trana, tranb;
	int m, n, k;
	double alpha = 1.0, beta = 1.0, StressForce[2];
	trana = 'N';
	tranb = 'N';
	int lda, ldb, ldc;
	NStress = Stress;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			StressForce[j] = 0;
		}
		
		for (int ixi = 0; ixi < 3; ixi++)
		{
			double xi = Rstg[ixi];
			double xih = H[ixi];
			for (int ieta = 0; ieta < 3; ieta++)
			{
				double eta = Rstg[ieta];
				double etah = H[ieta];
				double iStressForce[2];
				for (int j = 0; j < 2; j++)
				{
					iStressForce[j] = 0;
				}
				Shape(xi, eta);			
				B[0][0] = PDeri[0][i];
				B[0][1] = 0;
				B[1][0] = 0;
				B[1][1] = PDeri[1][i];
				B[2][0] = PDeri[1][i];
				B[2][1] = PDeri[0][i];
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						S[j][k] = 0;
					}
				}
				m = 2; n = 1; k = 3;
				lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
				dgemm(&trana, &tranb, &m, &n, &k, &alpha, B[0], &lda, NStress, &ldb, &beta, iStressForce, &ldc);
				for (int j = 0; j < 2; j++)
				{
					iStressForce[j] = iStressForce[j] * Det * xih *etah;
				}
				for (int j = 0; j < 2; j++)
				{
					StressForce[j] += iStressForce[j];
				}
			}
		}
		if (Dof[i*2] != 0)
		{
			InitialForce[Dof[i * 2] - 1] -= StressForce[0];
		}
		if (Dof[i * 2 + 1] != 0)
		{
			InitialForce[Dof[i * 2 + 1] - 1] -= StressForce[1];
		}
	}
	return 0;
}
int Quadrilateral::InitialStrain(double *InitialForce)
{
	double Elastic, Posion, Density;
	Material = Groups.Get(Group, "Material");
	Elastic = Mats.Get(Material, "Modulus");
	Posion = Mats.Get(Material, "Posion");
	Density = Mats.Get(Material, "Density");
	double D1 = Elastic*(1.000 - Posion) / ((1.00 + Posion)*(1 - 2.00*Posion));
	double D2 = Elastic*Posion / ((1.00 + Posion)*(1.00 - 2.00*Posion));
	double D3 = Elastic*0.50 / (1.00 + Posion);
	double D[3][3] = { D1, D2, 0, D2, D1, 0, 0, 0, D3 };
	double B[3][2], S[2][3] = { 0 };
	double *NStrain;
	char trana, tranb;
	int m, n, k;
	double alpha = 1.0, beta = 1.0, StrainForce[2];
	trana = 'N';
	tranb = 'N';
	int lda, ldb, ldc;
	NStrain = Strain;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			StrainForce[j] = 0;
		}

		for (int ixi = 0; ixi < 3; ixi++)
		{
			double xi = Rstg[ixi];
			double xih = H[ixi];
			for (int ieta = 0; ieta < 3; ieta++)
			{
				double eta = Rstg[ieta];
				double etah = H[ieta];
				double iStrainForce[2];
				for (int j = 0; j < 2; j++)
				{
					iStrainForce[j] = 0;
				}
				Shape(xi, eta);
				B[0][0] = PDeri[0][i];
				B[0][1] = 0;
				B[1][0] = 0;
				B[1][1] = PDeri[1][i];
				B[2][0] = PDeri[1][i];
				B[2][1] = PDeri[0][i];
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						S[j][k] = 0;
					}
				}
				m = 2; n = 3; k = 3;
				lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
				dgemm(&trana, &tranb, &m, &n, &k, &alpha, B[0], &lda, D[0], &ldb, &beta, S[0], &ldc);
				m = 2; n = 1; k = 3;
				lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
				dgemm(&trana, &tranb, &m, &n, &k, &alpha, S[0], &lda, NStrain, &ldb, &beta, iStrainForce, &ldc);
				for (int j = 0; j < 2; j++)
				{
					iStrainForce[j] = iStrainForce[j] * Det * xih *etah;
				}
				for (int j = 0; j < 2; j++)
				{
					StrainForce[j] += iStrainForce[j];
				}
			}
		}
		if (Dof[i * 2] != 0)
		{
			InitialForce[Dof[i * 2] - 1] += StrainForce[0];
		}
		if (Dof[i * 2 + 1] != 0)
		{
			InitialForce[Dof[i * 2 + 1] - 1] += StrainForce[1];
		}
	}
	return 0;
}
int Quadrilateral::ElementLoad(double *LoadMatrix)
{

	int m, n, k, lda, ldb, ldc;
	char trana='N', tranb='N';
	double alpha = 1.0, beta = 1.0;
	m = 8; n = 1; k = 8;
	lda = max(1, m); ldb = max(1, k); ldc = max(1, m);
	for (int i = 0; i < 8; i++)
	{
		Load[i] = 0;
	}
	dgemm(&trana, &tranb, &m, &n, &k, &alpha, Stiff[0], &lda, Result, &ldb, &beta, Load, &ldc);
	for (int i = 0; i < 8; i++)
	{
		if (Dof[i] != 0)
		{
			LoadMatrix[Dof[i] - 1] += Load[i];
		}		
	}
	for (int i = 0; i < 4; i++)
	{
		double NodeForce[2];
		NodeForce[0] = Load[2 * i];
		NodeForce[1] = Load[2 * i + 1];
		Nodes.PutResult(Node[i], NULL, "Time");
		Nodes.PutResult(Node[i], NodeForce, "Force");
	}
	return 0;
}
int Element::InitialStrain(double *InitialForce)
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].InitialStrain(InitialForce);
		}
	}
	return 0;
}
int Element::InitialStress(double *InitialForce)
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].InitialStress(InitialForce);
		}
	}
	return 0;
}
int Element::ElementLoad(double *LoadMatrix)
{
	int nelem = Type.size();
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		if (Type[ielem] == 4)
		{
			Quadrs[Index[ielem]].ElementLoad(LoadMatrix);
		}
	}
	return 0;
}
int Element::GetType(int ielem)
{
	return Type.at(ielem);
}
int Quadrilateral::GetGroup()
{
	return Group;
}
double Element::GetDensity(int ielem)
{
	if (Type[ielem] == 4)
	{
		int Group=Quadrs[Index[ielem]].GetGroup();
		int Material;
		double Density;
		Material = Groups.Get(Group, "Material");
		Density = Mats.Get(Material, "Density");
		return Density;
	}
	return 0;
}
