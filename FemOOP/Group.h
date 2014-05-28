#ifndef GROUP_H
#define GROUP_H
#pragma once
#include <string>
#include <vector>
#include <sstream>
using namespace std;
class Info
{
	int Index,Type, Nelem, Material;
	int *Elem;
public:
	Info(string text);
	int Get(string VarName);
	int *GetElem();
	int PutElem(int ielem,int ElemIndex);
};
class Group
{
	vector<Info> GList;
public:
	int ReadFile(ifstream &grp,int ngroup);
	int Get(int igroup,string VarName);
	int *GetElem(int igroup);
	int PutElem(int igroup,int ielem, int ElemIndex);
};


#endif