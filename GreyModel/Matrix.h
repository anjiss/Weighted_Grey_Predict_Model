#pragma once
#include <fstream>
using namespace std;

class Matrix
{
public:
	Matrix();
	Matrix(int n,int m);
	Matrix(int n, int m,double d);
	Matrix(int n, int m,double d[]);
	~Matrix();
	Matrix & operator +(Matrix &a);
	Matrix & operator -(Matrix &a);
	Matrix & operator *(Matrix &a);
	Matrix & operator *(double d);
	Matrix & Transpose();
	Matrix & Inverse();
	Matrix & Change(Matrix x);
	Matrix & operator= (const Matrix &a);
	void input();
	void disp();
	void outData(ofstream &myout);

	int row;//���������
	int col;//���������
	double **data;
};
