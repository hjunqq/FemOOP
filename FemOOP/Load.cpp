#include "Load.h"
#include <math.h>
extern Node Nodes;
Load *Loads;
extern Element Elems;
extern Group Groups;
const double Rstg[] = { -0.774596669241483, 0.0, 0.774596669241483 };
const double H[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
Load::Load()
{
}

int Load::ReadFile(ifstream &loa)
{
	string text;	
	int  iface = 0, ibody = 0, icent = 0;
	stringstream stream;
	
	getline(loa, text);
	getline(loa, text);
	stream << text;
	stream >> nloa >> load[0] >> load[1] >> load[2];
	stream.clear();
	Face = new PFace[load[0]];
	Body = new PBody[load[1]];
	Center = new PCenter[load[2]];
	for (int iloa = 0; iloa < nloa; iloa ++)
	{
		getline(loa, text);
		getline(loa, text);
		type.append(text+",");
		if (text == "face")
		{
			int nelem;
			getline(loa, text);
			getline(loa, text);
			Face[iface].SetValue(text);
			stream << text;
			stream >> nelem;
			stream.clear();
			getline(loa, text);
			for (int ielem = 0; ielem < nelem; ielem++)
			{				
				getline(loa, text);
				stream << text;
				Face[iface].SetFace(text, ielem);
				stream.clear();
			}
			iface++;
		}
		else if (text=="body")
		{
			getline(loa, text);
			getline(loa, text);
			Body[ibody].SetValue(text);
			getline(loa, text);
			getline(loa, text);
			Body[ibody].SetBody(text);
			ibody++;
		}
		else if (text=="cent")
		{
			getline(loa, text);
			getline(loa, text);
			Center[icent].SetValue(text);
			getline(loa, text);
			getline(loa, text);
			Center[icent].SetNode(text);
			icent++;
		}
	}
	return 0;
}
int PFace::SetValue(string text)
{
	stringstream stream;
	stream << text;
	stream >> Nelem >> Nnode >> Dir >> StartCoor >> EndCoor >> StartVal >> EndVal >> Curve;
	Dir--; Curve--;
	Elem = new int[Nelem];
	Node = new int*[Nelem];
	for (int ielem = 0; ielem < Nelem; ielem++)
	{
		Node[ielem] = new int[Nnode];
	}
	return 0;
}
int PFace::SetFace(string text, int iFace)
{
	stringstream stream;
	stream << text;
	for (int inode = 0; inode < Nnode; inode++)
	{
		stream >> Node[iFace][inode];
		Node[iFace][inode]--;
	}
	stream >> Elem[iFace];
	Elem[iFace]--;
	stream.clear();
	return 0;
}
int PBody::SetValue(string text)
{
	stringstream stream;
	stream << text;
	stream >> Dir >> Acceleration>>nBodyGroup;
	Dir--;
	BodyGroup = new int[nBodyGroup];
	return 0;
}
int PBody::SetBody(string text)
{
	stringstream stream;
	stream << text;
	for (int iGroup = 0; iGroup < nBodyGroup; iGroup++)
	{
		stream >> BodyGroup[iGroup];
		BodyGroup[iGroup]--;
	}
	return 0;
}
int PCenter::SetValue(string text)
{
	stringstream stream;
	stream << text;
	stream >> Nnode>>Dir>>Value;
	Dir--;
	Node = new int[Nnode];
	return 0;
}
int PCenter::SetNode(string text)
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
Load::~Load()
{
}
int Load::ApplyLoad(double *LoadMatrix, int *DegreeOfFreedom)
{
	for (int iface = 0; iface < load[0]; iface++)
	{
		Face[iface].ApplyLoad(LoadMatrix, DegreeOfFreedom);
	}
	for (int ibody = 0; ibody < load[1]; ibody++)
	{
		Body[ibody].ApplyLoad(LoadMatrix, DegreeOfFreedom);
	}
	for (int icent = 0; icent < load[2]; icent++)
	{
		Center[icent].ApplyLoad(LoadMatrix, DegreeOfFreedom);
	}
	return 0;
}
int PFace::ApplyLoad(double *LoadMatrix,int *DegreeOfFreedom)
{
	double **Coor, *Pval, Length, Normal[2], *Function, xi, xih, Det, *ElementLoad, PLoad;
	Coor = new double*[Nnode];
	Pval = new double[Nnode];
	Function = new double[Nnode];
	ElementLoad = new double[Nnode];

	for (int ielem = 0; ielem < Nelem; ielem++)
	{
		ElementLoad[0] = ElementLoad[1] = 0;
		for (int inode = 0; inode < Nnode; inode++)
		{
			Coor[inode] = Nodes.GetCoor(Node[ielem][inode]);
			Pval[inode] = (Coor[inode][Dir] - StartCoor) / (EndCoor - StartCoor)*(EndVal - StartVal) + StartVal;
		}
		Length = sqrt(pow((Coor[0][0] - Coor[1][0]), 2) + pow((Coor[0][1] - Coor[1][1]), 2));
		Normal[1] = (Coor[0][0] - Coor[1][0]) / Length;
		Normal[0] = (Coor[0][1] - Coor[1][1]) / Length;
		for (int ixi = 0; ixi < 3; ixi++)
		{
			xi = Rstg[ixi];
			xih = H[ixi];
			Function[0] = 0.5*(1 - xi);
			Function[1] = 0.5*(1 + xi);
			Det = 0.5*Length;
			PLoad = Function[0] * Pval[0] + Function[1] * Pval[1];
			ElementLoad[0] += Function[0] * PLoad * Det * xih;
			ElementLoad[1] += Function[1] * PLoad * Det * xih;
		}
		for (int inode = 0; inode < Nnode; inode++)
		{
			int NodeIndex;
			NodeIndex = Node[ielem][inode] ;
			if (DegreeOfFreedom[NodeIndex * 2] != 0)
			{
				LoadMatrix[DegreeOfFreedom[NodeIndex * 2]-1] += -ElementLoad[inode] * Normal[0];
			}
			if (DegreeOfFreedom[NodeIndex * 2 + 1] != 0)
			{
				LoadMatrix[DegreeOfFreedom[NodeIndex * 2 + 1]-1] += -ElementLoad[inode] * Normal[1];
			}
		}	
	}
	return 0;
}
int PBody::ApplyLoad(double *LoadMatrix, int *DegreeOfFreedom)
{
	double **Coor,**NLoad,*Function,Det,Density,*ELoad;
	int *Node, Etype, *Elem, GroupIndex,Nelem;
	for (int iGroup = 0; iGroup < nBodyGroup; iGroup++)
	{
		GroupIndex = BodyGroup[iGroup];
		Nelem = Groups.Get(GroupIndex, "Nelem");
		Elem = new int[Nelem];
		Elem = Groups.GetElem(GroupIndex);
		for (int ielem = 0; ielem < Nelem; ielem++)
		{
			Etype = Elems.GetType(Elem[ielem]);
			Density = Elems.GetDensity(ielem);
			if (Etype == 4)
			{
				Coor = new double*[4];
				Node = new int[4];
				Function = new double[4];
				ELoad = new double[2]();
				ELoad[Dir] = Density*Acceleration;
				NLoad = new double*[4];
				for (int inode = 0; inode < 4; inode++)
				{
					Node[inode] = Elems.GetNode(Elem[ielem], inode);
					Coor[inode] = Nodes.GetCoor(Node[inode]);
					NLoad[inode] = new double[2]();
				}
				for (int ixi = 0; ixi < 3; ixi++)
				{
					double xi = Rstg[ixi];
					double xih = H[ixi];
					for (int ieta = 0; ieta < 3; ieta++)
					{
						double eta = Rstg[ieta];
						double etah = H[ieta];
						Shape(xi, eta, Coor, Function, Det);
						for (int i = 0; i < 4; i++)
						{
							for (int idim = 0; idim < 2; idim++)
							{
								NLoad[i][idim] += Function[i] * ELoad[idim] * Det*xih*etah;
							}
						}
					}
				}
				for (int inode = 0; inode < 4; inode++)
				{
					for (int idim = 0; idim < 2; idim++)
					{
						if (DegreeOfFreedom[Node[inode] * 2 + idim] != 0)
						{
							LoadMatrix[DegreeOfFreedom[Node[inode] * 2 + idim] - 1] += NLoad[inode][idim];
						}
					}
				}
			}

		}
		delete[] Elem;
	}
	return 0;
}
int PCenter::ApplyLoad(double *LoadMatrix, int *DegreeOfFreedom)
{
	double NValue[2] = { 0 };
	NValue[Dir] = Value;
	for (int inode = 0; inode < Nnode; inode++)
	{
		for (int idim = 0; idim < 2; idim++)
		{
			if (DegreeOfFreedom[Node[inode] * 2 + idim] != 0)
			{
				LoadMatrix[DegreeOfFreedom[Node[inode] * 2 + idim] - 1] += NValue[idim];
			}
		}
	}
	return 0;
}
int PFace::Shape(double xi, double *Coor, double *Function, double *Deri)
{
	double DeriFunction[2],Xjac;
	Function[0] = 0.5*(1 - xi);
	Function[1] = 0.5*(1 + xi);
	DeriFunction[0] = -0.5;
	DeriFunction[1] = 0.5;
	Xjac = 0.0;
	for (int i = 0; i < 2; i++)
	{
		Xjac += DeriFunction[i] * Coor[i];
	}
	return 0;
}
int PBody::Shape(double xi,double eta, double **Coor, double *Function, double & Det)
{
	double PartialDeri[2][4], Xjac[2][2];
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
	return 0;
}