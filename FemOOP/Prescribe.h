#ifndef PRE_H
#define PRE_H
#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Element.h"
using namespace std;
class FixBoundary
{
	int Dir, Nnode;
	int *Node;
public:
	FixBoundary(string text);
	int AddVal(string text);
	int Get(string text);
	int Get(int index);
};
class Displacement
{
	int Dir, Nnode;
	double Value;
	int *Node;
public:
	Displacement(string text);
	int AddNode(string text);
	int ApplyDisp(double *LoadMatrix, int *DegreeOfFreedom, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement);
	int GetNnode();
	int GetDir();
	int GetNode(int inode);
};
class InteractBoundary
{
	int nNode, *Node,InteractDomain,*InteractNode; 
public:
	InteractBoundary(string text);
	int AddNode(string text);
	int AddInteractNode(string text);
	int GetNnode();
	int GetNode(int inode); 
	int GetInteractNode(int inode);
	int GetInteractDomain();
	int ApplyInterDisp(double *LoadMatrix, double  **StiffMatrix, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement, double *InteractResult, int *Interactnode);
};
class Prescribe
{
	int npre,ndim,nDisp,nInteract;
	vector<FixBoundary> Fix;
	vector<Displacement>Disp;
	vector<InteractBoundary>Interact;
	vector<int>Index;
	vector<string> type;
public:
	Prescribe();
	int ReadFile(ifstream &pre);
	int FixDof(int *DegreeOfFreedom, int n);
	int ApplyDisp(double *LoadMatrix, int *DegreeOfFreedom, int TotalDegreeOfFreedom, int TotalNode, int NodeDof,int nelem,double *InitialDisplacement );
	int ApplyInterDisp(double *LoadMatrix, double  **StiffMatrix, int TotalDegreeOfFreedom, int TotalNode, int NodeDof, int nelem, double *InitialDisplacement, double *InteractResult, int *Interactnode);
	int GetnInteract();
	int GetInteractProcess(int iInteract);
	int GetInteractnNode(int iInteract);
	int *GetInteractNode(int iInteract);
	~Prescribe();
};

#endif