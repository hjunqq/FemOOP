#include "Group.h"
#include "FEMOOP.h"
Group Groups;
Info::Info(string text)
{
	stringstream stream;
	stream<< text;
	stream >> Index >> Type >> Nelem >> Material;
	Index--;
	Material--;
	Elem = new int[Nelem];
}
int Info::Get(string VarName)
{
	if (VarName == "Type")
	{
		return Type;
	}
	else if (VarName == "Nelem")
	{
		return Nelem;
	}
	else if (VarName == "Material")
	{
		return Material;
	}
	return -1;
}
int *Info::GetElem()
{
	return Elem;
}
int Group::Get(int igroup,string VarName)
{
	return GList[igroup].Get(VarName);
}
int *Group::GetElem(int igroup)
{
	return GList[igroup].GetElem();
}
int Info::PutElem(int ielem, int ElemIndex)
{
	Elem[ielem] = ElemIndex;
	return 0;
}
int Group::PutElem(int iGroup, int ielem, int ElemIndex)
{
	GList[iGroup].PutElem(ielem, ElemIndex);
	return 0;
}
int Group::ReadFile(ifstream &grp,int ngroup)
{
	string text;
	stringstream stream;
	
	for (int igroup = 0; igroup < ngroup; igroup++)
	{
		getline(grp, text);
		getline(grp, text);
		GList.push_back(text);
	}
	return 0;
}
