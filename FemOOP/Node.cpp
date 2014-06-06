#include "Node.h"
#include "FEMOOP.h"
Node Nodes;


int Node::ReadFile(ifstream &cor,int ndim,int nnode)
{
	this->nnode = nnode;
	val = new double*[nnode];
	Disp = new double*[nnode];
	Strain = new double *[nnode];
	Stress = new double *[nnode];
	InternalForce = new double*[nnode];
	Nelem = new int[nnode]();
	for (int inode = 0; inode < nnode; inode++)
	{
		val[inode] = new double[ndim];
		Disp[inode] = new double[ndim]();
		Strain[inode] = new double[3]();
		Stress[inode] = new double[3]();
		InternalForce[inode] = new double[ndim]();
	}
	
	stringstream stream;
	string text;
	
	int idx;
	if (ndim == 2)
	{
		for (int inode = 0; inode < nnode; inode++)
		{
			getline(cor, text);
			stream << text << endl;
			stream >> idx >> val[inode][0] >> val[inode][1];
			idx--;
			stream.clear();
		}
	}
	else if (ndim == 3)
	{
		for (int inode = 0; inode < nnode; inode++)
		{
			getline(cor, text);
			stream << text << endl;
			stream >> idx >> val[inode][0] >> val[inode][1] >> val[inode][2];
			idx--;
			stream.clear();
		}
	}
	return 0;
}
double *Node::GetCoor(int inode)
{
	return val[inode]; 
}
double *Node::GetDisp(int inode)
{
	return Disp[inode];
}
double *Node::GetInternalForce(int inode)
{
	return InternalForce[inode];
}
double *Node::GetStrain(int inode)
{
	return Strain[inode];
}
double *Node::GetStress(int inode)
{
	return Stress[inode];
}
int Node::PutResult(int inode, double *Value,string ValType)
{
	if (ValType == "Disp")
	{
		for (int i = 0; i < 2; i++)
		{
			Disp[inode][i] = Value[i];
				//(Disp[inode][i] *( Nelem[inode] - 1) + Value[i]) / Nelem[inode];
		}
	}
	else if (ValType == "Time")
	{
		Nelem[inode]++;
	}
	else if (ValType == "Strain")
	{
		for (int i = 0; i < 3; i++)
		{
			this->Strain[inode][i]=
				(this->Strain[inode][i]*(Nelem[inode]-1) + Value[i])/Nelem[inode];

		}
	}
	else if (ValType == "Stress")
	{
		for (int i = 0; i < 3; i++)
		{
			this->Stress[inode][i] =
				(this->Stress[inode][i] * (Nelem[inode] - 1) + Value[i]) / Nelem[inode];

		}
	}
	else if(ValType == "Force")
	{
		for (int i = 0; i < 2; i++)
		{
			this->InternalForce[inode][i] =
				(this->InternalForce[inode][i] *(Nelem[inode] - 1) + Value[i])/ Nelem[inode];
		}
	}
	return 0;
}
int Node::PrintResult(ofstream &chk)
{
	for (int inode = 0; inode < nnode; inode++)
	{
		chk << setw(10) << inode
			<< setw(15) << Disp[inode][0]
			<< setw(15) << Disp[inode][1]<<endl;
	}
	return 0;
}
int Node::SetZero(string Type)
{
	if (Type == "All")
	{
		for (int inode = 0; inode < nnode; inode++)
		{
			for (int i = 0; i < 3; i++)
			{
				Stress[inode][i] = 0;
				Strain[inode][i] = 0;
				InternalForce[inode][i] = 0;
			}
			for (int i = 0; i < 2; i++)
			{
				Disp[inode][i] = 0;
			}
			Nelem[inode] = 0;
		}
	}
	else if (Type == "Times")
	{
		for (int inode = 0; inode < nnode; inode++)
		{
			Nelem[inode] = 0;
		}
	}
	
	return 0;
}