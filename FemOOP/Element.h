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

///四边形单元
class Quadrilateral 
{
	///单元序号
	int Index,		
		///单元节点个数
		Nnode=4,	
		///单元所属单元组
		Group,		
		///单元材料号
		Material,	
		///单元自由度个数
		Ndof,	
		///单元节点数组
		Node[4],	
		///单元自由度数组
		*Dof;		
		///单元刚度矩阵
	double Stiff[8][8],		
		///单元形函数
		Function[4],		
		///单元形函数的偏导数
		PartialDeri[2][4],	
		///单元节点坐标
		**Coor,				
		///Jacobi矩阵
		Xjac[2][2],			
		///Jacobi逆矩阵
		Rjac[2][2],			
		///B矩阵
		B[3][8],			
		///节点应变
		NodeStrain[4][3],	
		GaussStrain[4][3],
		///节点应力
		NodeStress[4][3],
		GaussStress[4][3],
		///Jacobi行列式的值
		Det,				
		///形函数对整体坐标的偏导数
		PDeri[2][4],		
		///节点位移
		Result[8] ,
		///单元应力
		Stress[3] ,
		///单元应变
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