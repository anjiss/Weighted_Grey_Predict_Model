// GreyModel.cpp : �������̨Ӧ�ó������ڵ㡣

#include "stdafx.h"
#include "Grey.h"
#include "Matrix.h"
#include <fstream>
using namespace std;

void main()
{
	ofstream myout("outdata.txt");
	Matrix X0,X1,X0_e,X0_temp;
	X0=GetTextToMatrix();
	X0_temp=GetTextToMatrix();
	myout<<"X0Ϊ��"<<endl;
	X0.outData(myout);
	myout<<endl;

	/*double *relation=new double[X0.col];
	myout<<"�ҹ�����:"<<endl;
	for (int j=0;j<X0.col;j++)
	{
		GreyRelation(X0_temp,j,relation);
		myout<<"��"<<j<<"������������������Ļҹ�����Ϊ";
		for (int i=0;i<X0.col;i++)
		{
			myout<<relation[i]<<"    ";
		}
		myout<<endl;
	}*/

	ResultX0(X0,X0.row+3,X0_e,3,myout);
//	ResultX0(X0_e,X0.row,X0_e,3,myout);
	double c,p;
	myout<<"Ԥ����Ϊ��"<<endl;
	X0_e.outData(myout);
	myout<<"�����"<<endl;
	for (int i=0;i<X0.col;i++)
	{
		GetError(c,p,X0,X0_e,i+1);
		myout<<c<<"    "<<endl;
	}


	myout.close();
}