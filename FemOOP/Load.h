#ifndef LOA_H
#define LOA_H
#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Node.h"
#include "Element.h"
#include "Group.h"

using namespace std;
class PFace
{
	double StartCoor, EndCoor, StartVal, EndVal;
	int Curve,Dir;
	int *Elem, **Node, Nelem, Nnode;
public:
	int SetValue(string stream);
	int SetFace(string stream, int iFace);
	int ApplyLoad(double *LoadMatrix, int *DegreeOfFreedom);
	int Shape(double xi, double *Coor, double *Function, double *Deri);
};
class PBody
{
	int Dir,nBodyGroup,*BodyGroup;
	double Acceleration;
public:
	int SetValue(string stream);
	int SetBody(string stream);
	int ApplyLoad(double *LoadMatrix, int *DegreeOfFreedom);
	int Shape(double xi,double eta, double **Coor, double *Function,double &Det);
};
class PCenter
{
	int Nnode,Dir,*Node;
	double Value;
public:
	int SetValue(string stream);
	int SetNode(string Stream);
	int ApplyLoad(double *LoadMatrix, int *DegreeOfFreedom);
};
class Load
{
	int nloa, load[3];
	string type;
	PFace *Face;
	PBody *Body;
	PCenter *Center;
public:
	Load();
	int ReadFile(ifstream &loa);
	int ApplyLoad(double *LoadMatrix,  int *DegreeOfFreedom);
	~Load();
};

#endif