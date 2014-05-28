#ifndef MAT_H
#define MAT_H
#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;
class Elastic
{
	double Modulus, Posion, Density;
public:
	Elastic(string text);
	double Get(string VarName);
};
class Material
{
	int mat;
	vector<string> type;
	vector<int> Index;
	vector<Elastic> LinearElastic;
public:
	Material();
	int ReadFile(string path, int nmat);
	double Get(int iMat, string VarName);
	~Material();
};

#endif