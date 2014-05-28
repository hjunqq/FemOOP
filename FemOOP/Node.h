#ifndef NODE_H
#define NODE_H
#pragma once
#include <string>

using namespace std;
class Node
{
	int nnode;
	double **val;
	double **Disp, **Strain, **Stress,**InternalForce;
	int *Nelem;
public:
	friend class Element;
	int ReadFile(ifstream &cor,int ndim,int nnode);
	double *GetCoor(int inode);
	double *GetStrain(int inode);
	double *GetDisp(int inode);
	double *GetStress(int inode);
	double *GetInternalForce(int inode);
	int PutResult(int inode, double *Value,string Valtype);
	int PrintResult(ofstream &chk);
	int SetZero(string Type);
};
#endif