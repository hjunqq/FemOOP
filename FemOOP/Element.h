#ifndef ELE_H
#define ELE_H
#pragma once
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "Node.h"
#include "Group.h"
#include "Material.h"

///�ı��ε�Ԫ
class Quadrilateral 
{
	///��Ԫ���
	int Index,		
		///��Ԫ�ڵ����
		Nnode=4,	
		///��Ԫ������Ԫ��
		Group,		
		///��Ԫ���Ϻ�
		Material,	
		///��Ԫ���ɶȸ���
		Ndof,	
		///��Ԫ�ڵ�����
		Node[4],	
		///��Ԫ���ɶ�����
		*Dof;		
		///��Ԫ�նȾ���
	double Stiff[8][8],		
		///��Ԫ�κ���
		Function[4],		
		///��Ԫ�κ�����ƫ����
		PartialDeri[2][4],	
		///��Ԫ�ڵ�����
		**Coor,				
		///Jacobi����
		Xjac[2][2],			
		///Jacobi�����
		Rjac[2][2],			
		///B����
		B[3][8],			
		///�ڵ�Ӧ��
		NodeStrain[4][3],	
		GaussStrain[4][3],
		///�ڵ�Ӧ��
		NodeStress[4][3],
		GaussStress[4][3],
		///Jacobi����ʽ��ֵ
		Det,				
		///�κ��������������ƫ����
		PDeri[2][4],		
		///�ڵ�λ��
		Result[8] ,
		///��ԪӦ��
		Stress[3] ,
		///��ԪӦ��
		Strain[3],
		Load[8] ;
public:
	Quadrilateral(string text);
	int EGetCoor();
	int Shape(double xi, double eta);
	int ElementStiff();
	double **GetStiff();
	int *GetDof();
	int GetGroup();
	int GetNode(int inode);
	int FillDof(int * DegreeOfFreedom);
	int GetResult(double *Result);
	int GetInitialDisplacement(double *InitialDisplacement);
	int GetInteractDisp(double *InteractDisp, int *InteractNode, int nInteractNode);
	int ElementStress();
	int PrintStress(ofstream & outdist);
	int PrintStress(ostream & outdist);
	int InitialStrain(double *InitialForce);
	int InitialStress(double *InitialForce);
	int ElementLoad(double *LoadMatrix);
};
class Element
{
	vector<int> Type;
	vector<int> Index;
	vector<Quadrilateral> Quadrs;
public:
	Element();
	int ReadFile(ifstream &ele, Group &Groups, int ngroup);
	int ElementStiff();
	double **GetStiff(int ielem);
	int *GetDof(int ielem);
	int GetNode(int ielem, int inode);
	int GetType(int ielem);
	double GetDensity(int ielem);
	int FillDof(int * DegreeOfFreedom);
	int GetResult(double *Result);
	int GetInitialDisplacement(double *InitialDisplacement);
	int GetInteractDisp(double *InteractDisp,int *InteractNode,int nInteractNode);
	int Ndof(int ielem);
	int ElementStress();
	int PrintStress(int ielem,ofstream & outdist);
	int PrintStress(int ielem, ostream & outdist);
	int InitialStrain(double *InitialForce);
	int InitialStress(double *InitialForce);
	int ElementLoad(double *LoadMatrix);
};

#endif