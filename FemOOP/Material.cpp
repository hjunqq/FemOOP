#include "Material.h"

Material Mats;
Material::Material()
{
}
Elastic::Elastic(string text)
{
	stringstream stream;
	stream << text;
	stream >> Modulus >> Posion >> Density;
}
int Material::ReadFile(string path, int nmat)
{
	string text;
	stringstream stream;
	ifstream mat(path + ".mat");
	if (mat.fail())
	{
		cout << "The material file does not exist" << endl;
		cout << "Exit program" << endl;
		return -1;
	}
	for (int imat = 0; imat < nmat; imat++)
	{
		getline(mat, text);
		stream << text;
		type.push_back(text);
		if (text == "Elastic")
		{
			getline(mat,text);
			LinearElastic.push_back(Elastic(text));
		}
	}
	return 0;
}
double Elastic::Get(string VarName)
{
	if (VarName == "Modulus")
	{
		return Modulus;
	}
	else if (VarName == "Posion")
	{
		return Posion;
	}
	else if (VarName == "Density")
	{
		return Density;
	}
	return -1;
}
double Material::Get(int iMat, string VarName)
{
	if (type.at(iMat) == "Elastic")
		return LinearElastic[iMat].Get(VarName);
	return -1;
}
Material::~Material()
{
}
