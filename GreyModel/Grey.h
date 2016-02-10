#pragma once

#include "StdAfx.h"
#include "string.h"
#include "Matrix.h"
#include "math.h"
//#include "Grey.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
using namespace std;


Matrix & GetX1(Matrix & X0)//累加得到X1
{
	int n=X0.col;
	int m=X0.row;
	Matrix *this_X1=new Matrix(m,n);
	for (int j=0;j<n;j++)
	{
		this_X1->data[0][j]=X0.data[0][j];
	}
	for (int i=1;i<m;i++)
	{
		for (int j=0;j<n;j++)
		{
			this_X1->data[i][j]=X0.data[i][j]+this_X1->data[i-1][j];
		}
	}
	return *this_X1;
}

Matrix & GetL(Matrix & X1)//得到MGM建模过程中的L
{
	int n=X1.col;
	int m=X1.row;
	Matrix *this_L=new Matrix(m-1,n+1);
	for (int i=0;i<m-1;i++)
	{
		for (int j=0;j<n;j++)
		{
			this_L->data[i][j]=(X1.data[i][j]+X1.data[i+1][j])*0.5;
		}
		this_L->data[i][n]=1;
	}

	return *this_L;
}

Matrix & GetTextToMatrix()//将文件夹下的数据调入到程序中并以矩阵形式返回
{	
	Matrix *matx;

	FILE *file=fopen("indata.txt","r");
	int x,y;


	if (file==NULL)
	{
		cout<<"源文件打开失败！"<<endl;
	}
	else
	{
		fscanf(file,"%d,%d,\n",&x,&y);
		matx=new Matrix(x,y);
		for (int i=0;i<x;i++)
		{
			for (int j=0;j<y;j++)
			{
				fscanf(file,"%lf,",&(matx->data[i][j]));
			}
		}
	}
	
	return *matx;
}


Matrix & GetY(Matrix & X0)//得到MGM建模过程中的Y
{
	Matrix *this_Y=new Matrix(X0.row-1,X0.col);
	
	for (int j=0;j<X0.col;j++)
	{
		for(int i=0;i<X0.row-1;i++)
			{
				this_Y->data[i][j]=X0.data[i+1][j];
			}
	}
	
	return *this_Y;
}

void GetAB(Matrix &A,Matrix &B,Matrix &L,Matrix &Y)//利用L、Y算得建模过程中的A、B
{
	Matrix a,b;
	b=(L.Transpose()*L).Inverse();
	a=L.Transpose()*Y;
	a=b*a;

	Matrix m(a.col,a.row-1);
	Matrix n(a.col,1);
	for (int i=0;i<m.row;i++)
	{
		for (int j=0;j<m.col;j++)
		{
			m.data[i][j]=a.data[j][i];
		}
	}
	for (int i=0;i<n.row;i++)
	{
		n.data[i][0]=a.data[a.row-1][i];
	}
	A=m;
	B=n;
}
void Gete(Matrix &A,int t,Matrix &e)//利用参数A、t计算得到e
{
	Matrix temp(A.col,A.col,1);
	Matrix a=A*t;
	temp=temp+a;
	for (int i=2;i<100;i++)
	{
		a=a*A;
		a=a*(t/double(i));
		temp=temp+a;
	}
	e=temp;
}
void GreyRelation(Matrix &X0,int a,double *re,double f=0.5)//改正
{
	double **init=new double *[X0.col];
	double **minus=new double *[X0.col];
	double minus_min=0,minus_max=0;
	for (int i=0;i<X0.col;i++)
	{
		init[i]=new double[X0.row];
		minus[i]=new double[X0.row];
	}
	for (int i=0;i<X0.col;i++)
	{
		for (int j=0;j<X0.row;j++)
		{
			init[i][j]=X0.data[j][i]/X0.data[0][i];
		}
	}
	for (int i=0;i<X0.col;i++)
	{
		for (int j=0;j<X0.row;j++)
		{
			minus[i][j]=fabs(init[i][j]-init[a][j]);
			if (minus[i][j]>minus_max)
			{
				minus_max=minus[i][j];
			}
		}
	}
	double temp;
	for (int i=0;i<X0.col;i++)
	{
		temp=0;
		for (int j=0;j<X0.row;j++)
		{
			temp+=(minus_max*f)/(minus[i][j]+minus_max*f);
		}
		re[i]=temp/X0.row;
	}
	for(int i=0;i<X0.col;i++)
	{
		delete []minus[i];
		minus[i]=NULL;
	}
	delete [] minus;
	minus=NULL;
	for(int i=0;i<X0.col;i++)
	{
		delete []init[i];
		init[i]=NULL;
	}
	delete [] init;
	init=NULL;
}
Matrix & ajs(Matrix &X0,double *p)
{
	Matrix *X0_r=new Matrix(X0.row,X0.col);
	double c;
	for (int i=0;i<X0.row;i++)
	{
		X0_r->data[0][i]=X0.data[0][i];
	}
	for (int i=1;i<X0.row;i++)
	{
		for (int j=0;j<X0.col;j++)
		{
			c=X0.data[i][j]-X0.data[i-1][j];
			X0_r->data[i][j]=X0.data[i-1][j]+c*p[j];
		}
	}
	return *X0_r;
}

void ResultX0(Matrix &X0,int t,Matrix &X_re,int choise,ofstream &outdata)//给出一个已知矩阵X0，利用GMG模型算出预测矩阵，结果存放到X_re中
{
	//cout<<endl<<"X0"<<endl;
	//X0.disp();
//	Matrix P(X0.col+1,X0.col+1,1);
	Matrix X1,L,Y,A,B;
	Matrix X0_p;
	X0_p = X0;
	/*double *relation=new double[X0.col];
	GreyRelation(X0,choise,relation);
	X0_p=ajs(X0,relation);*/
	outdata<<"加权后的X0"<<endl;
	
	X0_p.outData(outdata);
	X1=GetX1(X0_p);
	L=GetL(X1);
	Y=GetY(X0_p);
	GetAB(A,B,L,Y);
	X0.disp();
	Matrix e;
	Matrix res1(t,X0.col);
	Matrix res0(t,X0.col);
	Matrix I(A.col,A.col,1);
	Matrix x(1,X0.col);//临时存储数据的单独一行
	Matrix x1(1,X0.col);
	for (int i=0;i<X1.col;i++)
	{
		x1.data[0][i]=X1.data[0][i];
	}
	for (int i=0;i<t;i++)
	{
		Gete(A,i,e);
		Matrix u1,u2;
		u1=e*x1.Transpose();
		u2=A.Inverse()*(e-I);
		u2=u2*B;
		x=u1+u2;
		for (int j=0;j<X1.col;j++)
		{
			res1.data[i][j]=x.Transpose().data[0][j];
		}
	}
	for (int j=0;j<res1.col;j++)
	{
		res0.data[0][j]=res1.data[0][j];
	}
	for (int i=1;i<t;i++)
	{
		for (int j=0;j<res1.col;j++)
		{
			res0.data[i][j]=res1.data[i][j]-res1.data[i-1][j];
		}
	}
	X_re=res0;
//	delete []relation;
//	relation=NULL;
}
void GetError(double &c,double &p,Matrix &x,Matrix &x_e,int t)//精度评定
{
	double *e=new double[x.row];
	double e_ave(0),x_ave(0);

	for (int i=0;i<x.row;i++)
	{
		e[i]=x.data[i][t-1]-x_e.data[i][t-1];
	}
	for (int i=0;i<x.row;i++)
	{
		e_ave+=e[i];
		x_ave+=x_e.data[i][t-1];
	}
	e_ave=e_ave/x.row;
	x_ave=x_ave/x.row;
	double S1(0),S2(0);
	for (int i=0;i<x.row;i++)//计算S1S2
	{
		S1+=(x.data[i][t-1]-x_ave)*(x.data[i][t-1]-x_ave);
		S2+=(e[i]-e_ave)*(e[i]-e_ave);
	}
	S1=S1/x.row;
	S2=S2/x.row;
	double C;
	C=sqrt(S2/S1);
	c=C;
	delete []e;
}




//void ResultX0WithP(Matrix &X0,int t,Matrix &X_re)
//{
//	Matrix res0(t,X0.col);
//	for (int ajs=0;ajs<X0.col;ajs++)
//	{
//		Matrix X1,L,Y,A,B,X0_temp,P1(X0.col+1,X0.col+1),P2(X0.col,X0.col),temp(X0.col,X0.col);
//		X0_temp=X0;
//		double *p=new double[X0.col];
//		GreyRelation(X0_temp,ajs,p);
//		
//		for (int i=0;i<X0.col;i++)
//		{
//			P1.data[i][i]=p[i];
//			P2.data[i][i]=p[i];
//		}
//		P1.data[X0.col][X0.col]=1;
//
//		X1=GetX1(&X0_temp);
//		L=GetL(&X1);
//		Y=GetY(&X0_temp);
//		GetAB(A,B,L,Y,P1);
//		Matrix e;
//		Matrix res1(t,X0.col);
//		Matrix I(A.col,A.col,1);
//		Matrix x(1,X0.col);//临时存储数据的单独一行
//		Matrix x1(1,X0.col);
//		for (int i=0;i<X1.col;i++)
//		{
//			x1.data[0][i]=X1.data[0][i];
//		}
//		for (int i=0;i<t;i++)
//		{
//			temp=P2;
//			Gete(A,i,e,temp);
//			Matrix u1,u2;
//			u1=e*x1.Transpose();
//			u2=A.Inverse();
//			u2=u2*P2;
//			u2=u2*(e-I);
//			u2=u2*B;
//			x=u1+u2;
//			for (int j=0;j<X1.col;j++)
//			{
//				res1.data[i][j]=x.Transpose().data[0][j];
//			}
//		}
//		res0.data[0][ajs]=res1.data[0][ajs];
//		for (int i=1;i<t;i++)
//		{
//			res0.data[i][ajs]=res1.data[i][ajs]-res1.data[i-1][ajs];
//		}
//
//	}
//	X_re=res0;
//}

Matrix & GetLWithP(Matrix & X1,Matrix &P)//得到MGM建模过程中的L
{
	int n=X1.col;
	int m=X1.row;
	Matrix *this_L=new Matrix(m-1,n+1);
	for (int i=0;i<m-1;i++)
	{
		for (int j=0;j<n;j++)
		{
			this_L->data[i][j]=(X1.data[i][j]+X1.data[i+1][j])*0.5;
		}
		this_L->data[i][n]=1;
	}

	return *this_L;
}

Matrix & GetP(Matrix &X0)
{
	Matrix *P=new Matrix(X0.col,X0.col);
	double *relation=new double[X0.col];
	for (int i=0;i<X0.col;i++)
	{
		GreyRelation(X0,i,relation);
		for (int j=0;j<X0.col;j++)
		{
			P->data[i][j]=relation[j];
		}
	}
	return *P;
}

