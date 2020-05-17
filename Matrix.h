#ifndef MATRIX
#define MATRIX

#include <iostream>
#include <ctime>
#include <vector>
using namespace std;

class Matrix {

private:
	double** matrix;
	int rows, cols;
	void _mStep(int row = 0, int col = 0);
	void _mRref(int r = 0, int c = 0);
	void _iStep(Matrix& inverse, int row = 0, int col = 0);
	void _iRref(Matrix& inverse, int r = 0, int c = 0);
	double _det(Matrix& m);

public:
	Matrix(int dim = 1, const double& init = double());

	Matrix(int dim, double* aij);
	Matrix(vector<int> v);

	Matrix(int rows, int cols, const double& init = double());
	Matrix(int rows, int cols, double** matrix);
	Matrix(const Matrix& m);
	bool isNull();
	bool isSquare();
	bool isDiagonal();
	bool isIdentity();
	bool isStep();
	bool isRref();
	bool isInvertible();
	
	bool isZeros();

	int size();

	void mSwap(int x, int y);
	void mMult(int x, double k);
	void mComb(int x, int y, double k = double(1));

	void print();

	Matrix swap(int x, int y);
	Matrix mult(int x, double k);
	Matrix comb(int x, int y, double k = double(1));
	Matrix transpost();

	void updateDiagonalMatrix();

	void mStep();
	Matrix step();
	void mRref();
	Matrix rref();
	Matrix mInvert();
	Matrix invert();
	int rang();
	Matrix sub(int row, int col);
	double det();
	Matrix& operator=(const Matrix& m);
	Matrix operator-();
	Matrix& operator++();
	Matrix& operator--();
	Matrix operator++(int n);
	Matrix operator--(int n);
	Matrix operator+(const Matrix& matrix);
	Matrix operator-(const Matrix& matrix);
	Matrix operator*(const Matrix& matrix);
	Matrix operator*(double k);
	Matrix operator/(double k);
	Matrix operator|(const Matrix& matrix);
	Matrix operator+=(const Matrix& matrix);
	Matrix operator-=(const Matrix& matrix);
	Matrix operator*=(double k);
	Matrix operator/=(double k);
	bool operator==(const Matrix& matrix);
	bool operator!=(const Matrix& matrix);
	double* operator[] (int i);
	friend Matrix operator*(double k, const Matrix& m);
	friend bool operator<(const Matrix& matrix, const Matrix& matrix2);
	friend ostream& operator<< (ostream& out, const Matrix& matrix);
	friend istream& operator>> (istream& in, Matrix& matrix);


};

Matrix operator*(double k, const Matrix& m);
Matrix diagonal(int n, double k);
Matrix identity(int n);

#endif
