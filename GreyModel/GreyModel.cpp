// GreyModel.cpp : 定义控制台应用程序的入口点。

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
	myout<<"X0为："<<endl;
	X0.outData(myout);
	myout<<endl;

	/*double *relation=new double[X0.col];
	myout<<"灰关联度:"<<endl;
	for (int j=0;j<X0.col;j++)
	{
		GreyRelation(X0_temp,j,relation);
		myout<<"第"<<j<<"个变量相对其他变量的灰关联度为";
		for (int i=0;i<X0.col;i++)
		{
			myout<<relation[i]<<"    ";
		}
		myout<<endl;
	}*/

	ResultX0(X0,X0.row+3,X0_e,3,myout);
//	ResultX0(X0_e,X0.row,X0_e,3,myout);
	double c,p;
	myout<<"预测结果为："<<endl;
	X0_e.outData(myout);
	myout<<"后验差"<<endl;
	for (int i=0;i<X0.col;i++)
	{
		GetError(c,p,X0,X0_e,i+1);
		myout<<c<<"    "<<endl;
	}


	myout.close();
}