#include "StdAfx.h"
#include "Matrix.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

Matrix::Matrix()
{
	row=0;
	col=0;
	data=NULL;
}
Matrix::Matrix(int n, int m)//���캯�����������N*M ����ʼ��Ϊ0����
{
	row=n;
	col=m;
	data=new double *[n];
	for(int i=0;i<n;i++)
	{
		data[i]=new double[m];
	}
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			data[i][j]=0.0;
}
Matrix::Matrix(int n, int m,double d)
{
	row=n;
	col=m;
	data=new double *[n];
	for(int i=0;i<n;i++)
	{
		data[i]=new double[m];
	}
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		{
			if (i!=j)
				data[i][j]=0.0;
			else
				data[i][j]=d;
		}
}

Matrix::Matrix(int n, int m,double d[])
{
	row=n;
	col=m;
	data=new double *[n];
	for(int i=0;i<n;i++)
	{
		data[i]=new double[m];
	}
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		{
			if (i!=j)
				data[i][j]=0.0;
			else
				data[i][j]=d[i];
		}

}

Matrix::~Matrix()
{
	for(int i=0;i<row;i++)
	{
		delete []data[i];
		data[i]=NULL;
	}
	delete [] data;
	data=NULL;
}

Matrix & Matrix::operator +(Matrix &a)
{
	if(row!=a.row||col!=a.col)
	{
		cout<<"\n\t\t��ͬ�о��󣬲������!\n";
		Matrix *c=new Matrix();
		return *c;
	}
	Matrix *b=new Matrix(a.row,a.col);
	for(int i=0;i<a.row;i++)
	{
		for(int j=0;j<a.col;j++)
		{
			b->data[i][j]=data[i][j]+a.data[i][j];
		}
	}
	return *b;
}

Matrix & Matrix::operator -(Matrix &a)
{
	if(row!=a.row||col!=a.col)
	{
		cout<<"\n\t\t��ͬ�о��󣬲������!\n";
		Matrix *c=new Matrix();
		return *c;
	}
	Matrix *b=new Matrix(a.row,a.col);
	for(int i=0;i<a.row;i++)
	{
		for(int j=0;j<a.col;j++)
		{
			b->data[i][j]=data[i][j]-a.data[i][j];
		}
	}
	return *b;
}

Matrix & Matrix::operator *(Matrix &a)
{
	if(col!=a.row)
	{
		cout<<"\n\t\t������������Ҫ�󣬲������!\n";
		Matrix *c=new Matrix();
		return *c;
	}
	Matrix *b=new Matrix(row,a.col);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<a.col;j++)
		{
//			b->data[i][j]=data[i][j]-a.data[i][j];
			double temp=0;
			for(int k=0;k<col;k++)
			{
				temp+=(data[i][k]*a.data[k][j]);
			}
			b->data[i][j]=temp;
		}
	}
	return *b;
}

Matrix & Matrix::operator *(double d)
{
	Matrix *b=new Matrix(row,col);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			b->data[i][j]=data[i][j]*d;
		}
	}
	return *b;

}
Matrix & Change(Matrix x)
{
	Matrix *b=new Matrix(x.row,x.col);
	for (int i=0;i<x.row;i++)
	{
		for (int j=0;j<x.col;j++)
		{
			b->data[i][j]=x.data[i][j];
		}
	}
	return *b;
}

Matrix & Matrix::Transpose()//ת��
{
	Matrix *b=new Matrix(col,row);
	for(int i=0;i<col;i++)
		for(int j=0;j<row;j++)
			b->data[i][j]=data[j][i];
	return *b;
}

Matrix & Matrix::Inverse()//����
{
	if(col!=row)//��������Ƿ���
	{
		cout<<"\n\t\t�þ�������������ȣ�������������!\n";
		Matrix *c=new Matrix();
		return *c;
	}
	Matrix *b=new Matrix(col,row);
	for(int k=0;k<col;k++)
	{
		Matrix i_M(col,col);
		i_M=*this;
		Matrix i_L(col,1);
		i_L.data[k][0]=-1;
//		cout<<"IM��IL�ĳ�ʼֵ"<<endl;
//		i_M.disp();
//		i_L.disp();
		for(int i=0;i<row-1;i++)
		{
			Matrix i_R(col,col,1);//
			for(int j=i+1;j<col;j++)
			{
				i_R.data[j][i]=-i_M.data[j][i]/i_M.data[i][i];
			}
//			cout<<"IRѭ��"<<i<<"�κ�Ľ��"<<endl;
//			i_R.disp();
			i_M=i_R*i_M;
			i_L=i_R*i_L;
//			cout<<"��ʱIM��ֵΪ"<<endl;
//			i_M.disp();
		}

		b->data[row-1][k]=-i_L.data[row-1][0]/i_M.data[row-1][row-1];
		for(int i=row-1;i>-1;i--)
		{
			double temp=0;
			for(int j=row-1;j>i;j--)
			{
				temp+=(b->data[j][k]*i_M.data[i][j]);
			}
			b->data[i][k]=(-i_L.data[i][0]-temp)/i_M.data[i][i];
		}
	}
	return *b;
}

Matrix & Matrix::operator = (const Matrix &a)
{
	if(data!=NULL)
	{
		for(int i=0;i<row;i++)
		{
			delete [] data[i];
			data[i]=NULL;
		}
		delete [] data;
		data=NULL;
	}
	row=a.row;
	col=a.col;
	data=new double *[row];
	for(int i=0;i<row;i++)
	{
		data[i]=new double[col];
	}
	for(int i=0;i<a.row;i++)
	{
		for(int j=0;j<a.col;j++)
		{
			data[i][j]=a.data[i][j];
		}
	}
	return *this;
}
void Matrix::input()
{
	cout<<"����������е����ݣ�";
	for(int i=0;i<row;i++)
	{
		cout<<"\n�������"<<i+1<<"�����ݣ�\n";
		for(int j=0;j<col;j++)
		{
			cout<<"data["<<i<<"]["<<j<<"]=";
			cin>>data[i][j];
		}
	}
}
void Matrix::disp()
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
			cout<<setprecision(12)<<data[i][j]<<'\t';
		cout<<endl;
	}
}

void Matrix::outData(ofstream &myout)
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
			myout<<data[i][j]<<'\t';
		myout<<endl;
	}
}