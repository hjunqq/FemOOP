#include "FEMOOP.h"

//////////////////////////////////////////////////////////////////////////  
///     COPYRIGHT NOTICE  
///     Copyright (c) 2014, �Ӻ���ѧ 
///     All rights reserved.  
///  
/// @file main.cpp
/// @brief  ����ƽ���ı��ε�Ԫ���
///  
///�����ļ�ʵ�ֵĹ��ܵ�������  
///  
/// @version 0.0
/// @author ��۾�
/// @date  2014��4��18��
///  
///  
/// 
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	FEMOOP Prob;
	string probn, text;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int n, myid, numprocs, i,namelen,mpi_comm;
	MPI::Init(argc, argv);
	numprocs = MPI::COMM_WORLD.Get_size();
	myid = MPI::COMM_WORLD.Get_rank();


	Prob.SendID(numprocs, myid);

	MPI::Get_processor_name(processor_name, namelen);
	cout << "Process " << myid << " of " << numprocs << " is on " <<
		processor_name << endl;
	ifstream inp("1.inp");
	if (inp.fail())
	{
		cout << "The input file does not exist" << endl;
		cout << "Exit program" << endl;
		return -1;
	}
	getline(inp, text);
	getline(inp, text);
	getline(inp, probn);
	if (myid == 0)
	{
		Prob.SetProbn(text);
	}
	if (myid == 1)
	{
		Prob.SetProbn(probn);
	}
	inp.close();

	Prob.StepCycle();

	MPI::Finalize();
}