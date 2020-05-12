// ErmitConstant2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Matrix.h"
#include <set>
#include <algorithm>

using namespace std;

const int N = 4;
const int MAXN = N * (N - 1) / 2;
double a[MAXN]; //used to generate matrix
double zeroArr[MAXN];
vector<Matrix> generatedVectors; // that x^T F x < 1
int genVect[N];
int v[N];

//todo create recursive generate vector for > 0 case
double history; //for higherOneCase

void generateMatrix(int n, int m, int ind, vector<Matrix>* globMatr) {
	if (!n && !m)
	{
		globMatr->push_back(Matrix(N, a));
		return;
	}

	if (n)
	{
		a[ind] = -0.5;
		generateMatrix(n - 1, m, ind + 1, globMatr);
	}

	if (m)
	{
		a[ind] = 0.5;
		generateMatrix(n, m - 1, ind + 1, globMatr);
	}
}

Matrix findMinMatrix(vector<Matrix>* matrixs) {
	Matrix min;
	for (Matrix m : *matrixs) {
		if (min.isNull()) {
			min = m;
		}
		else {
			if ((m.det() < min.det()) && (m.det() > 0)) {
				min = m;
			}
		}
	}
	return min;
}

double matMult(Matrix vect, Matrix matr) {
	Matrix m1 = (vect.transpost() * matr) * vect;
	return m1.det();
}

Matrix calculateNewMatrix(Matrix f, Matrix h, Matrix vector) {

	double fx = matMult(vector, f);
	double hx = matMult(vector, h);
	double c1 = 1 / (fx - hx);
	double c2 = 1 - hx;
	double c3 = fx - 1;
	Matrix tmp = c1 * ((c2 * f) + (c3 * h));
	return tmp;
}

/*
void recursiveGenVector(vector<int> v, int dim, vector<Matrix>* vm) {
	int k = N - dim;
	if (dim == 1) {
		for (int i = (-1) * k; i < 0; i++) {
			v.push_back(i);
			vm->push_back(Matrix(v));
			v.pop_back();
		}
	}
	else {
		for (int i = (-1) * k; i < 1 * k; i++) {
			v.push_back(i);
			recursiveGenVector(v, dim - 1, vm);
			v.pop_back();
		}
	}

}
*/
bool isAllZerosVector(vector<int> v) {
	for (int i : v) {
		if (i != 0) return false;
	}
	return true;
}

/*
vector<Matrix> generateVector() {
	if (N >= 2) {
		vector<Matrix> vm;
		vector<int> v;
		for (int i = (-1) * realxn; i <= realxn; i++) {
			xn = i;
			std::cout << "genVect[N] = " << genVect[N-1] << endl;
			for (int j = 0; j < sizeof(genVect) / sizeof(genVect[0]); j++) {
				int tmp = genVect[j];
				std::cout << "v[j] = " << tmp << endl;
				v.push_back(tmp);
			}
			if (!isAllZerosVector(v)) {
				vm.push_back(Matrix(v));
			}
			v.clear();
		}
		return vm;
	}
	else
		cerr << "dimention must be >= 2" << endl;
}

*/

/*
vector<Matrix> generateVector() {
	if (N >= 2) {
		vector<Matrix> vm;
		vector<int> v;
		for (int i = (-1) * genVect[0]; i < genVect[0]; i++) {
			v.push_back(i);
			recursiveGenVector(v, N - 1, &vm);
			v.pop_back();
		}
		return vm;
	}
	else
		cerr << "dimention must be >= 2" << endl;
}*/

bool recGenerateVectors(Matrix d, Matrix p, int n);

bool isAllZeros() {
	for (int j = 0; j < sizeof(genVect) / sizeof(genVect[0]); j++) {
		if (genVect[j] != 0) return false;
	}
	return true;
}

void generateVectorsFor(Matrix f) {
	generatedVectors.clear();
	Matrix e = Matrix(N, zeroArr);
	bool isNeededGenVect = recGenerateVectors(f, e, 0);
	//create vertor/vectors from int to Matrix
	if (isNeededGenVect) {
		std::cout << "Is needed to generate vectors" << endl;
		if (!isAllZeros()) {
			//generatedVectors = generateVector();
		}
		else {
			std::cout << "All zeros vector" << endl;
		}
	}
	else {
		vector<int> vect;
		for (int i = 0; i < N; i++) {
			vect.push_back(genVect[i]);
		}
		generatedVectors.push_back(Matrix(vect));
	}
	
}

double sumOfDiagonalElements(Matrix f, int n) {
	double counter = 0;
	for (int i = 0; i < n; i++) {
		counter += f[i][i];
	}
	return counter;
}

void reverseRecGenVectorForZeroCase(Matrix d, Matrix p, int n) {
	if (n > 0) {
		//востановить
		for (int i = n; i < d.size(); i++) {
			double coof = p[n-1][i] * (-1);
			d = d.comb(i, n-1, coof);
			p = p.comb(i, n-1, coof);
		}
		d.updateDiagonalMatrix();
		//посчитать элемент вектора n-1
		double coof = d[n - 1][n - 1];
		double xi = 0;
		for (int i = n; i < d.size(); i++) {
			xi -= d[n-1][i] / coof;
		}
		genVect[n - 1] = (int) xi;
		//рекурсивный вызов
		reverseRecGenVectorForZeroCase(d, p, n - 1);
	}
}

void higherZeroCaseRecFunc(Matrix d, Matrix p, int n, vector<int> vReverce, vector<Matrix>* vm) {
	if (n > 0) {
		//востановить
		cout << "d.size - n-1 = " << d.size() - n - 1 << endl;

		//bad?
		for (int i = n+1; i < d.size(); i++) {
			double coof = p[n][i] * (-1);
			d = d.comb(i, n, coof);
			p = p.comb(i, n, coof);
		}
		//

		d.updateDiagonalMatrix();
		//посчитать элемент вектора n-1

		//
		double a = 0; //bad
		for (int i = n+1; i < d.size(); i++) {
			a += d[n][i] * v[i];
		}
		//
		double b = 1 - history;
		double c = p[n][n]; //bad?
		double xi = sqrt(abs(b / c)) - a;
		cout << "history = " << history << endl;
		cout << "xi = " << xi << endl;
		genVect[n] = abs((int)xi);
		double coof = d[n][n];
		history += coof * (pow(v[n] + a, 2));
		//рекурсивный вызов
		for (int k = (-1) * genVect[n]; k <= genVect[n]; k++) {
			vReverce.push_back(k);
			v[n] = k;
			cout << "if (n > 0) k = " << k << ", v[" << n << "] = " << v[n] << endl;
			higherZeroCaseRecFunc(d, p, n - 1, vReverce, vm);
			vReverce.pop_back();
		}
	}
	else {
		for (int k = (-1) * genVect[n]; k <= genVect[n]; k++) {
			vReverce.push_back(k);
			v[n] = k;
			cout << "else k = " << k << ", v[" << n << "] = " << v[n] << endl;
			std::reverse(vReverce.begin(), vReverce.end());
			cout << endl;
			vm->push_back(Matrix(vReverce));
			vReverce.pop_back();
		}
	}
}

void reverceRecGenVectorForHigherZeroCase(Matrix d, Matrix p, int n) {
	//TODO: Create a recursive generate vector
	if (n > 0) {

		vector<Matrix> vm;
		vector<int> vReverce;
		for (int k = (-1) * genVect[n]; k <= genVect[n]; k++) {
			vReverce.push_back(k);
			v[n] = k;
			cout << "vm k = " << k << ", v[" << n << "] = " << v[n] << endl;
			higherZeroCaseRecFunc(d, p, n-1, vReverce, &vm);
			vReverce.pop_back();

		}
		generatedVectors = vm;
	}
}

bool recGenerateVectors(Matrix d, Matrix p, int n) {
	if (d[n][n] > 0) {
		if (n < d.size()-1) {
			//преобразования....
			for (int i = n+1; i < d.size(); i++) {
				double coof = (d[n][i] / d[n][n]) * (-1);
				d = d.comb(i, n, coof);
				p = p.comb(i, n, coof);
			}
			//
			d.updateDiagonalMatrix();
			return recGenerateVectors(d,p,n+1);
		}
		else {
			genVect[n] = abs(1 / (double)d[n][n]);
			v[n] = genVect[n];
			cout << "1 / (double)d[n][n] = " << 1 / (double)d[n][n] << " v["<< n << "] = " << v[n] << endl;
			history = pow(genVect[n], 2)*d[n][n];
			reverceRecGenVectorForHigherZeroCase(d, p, n);
			return true;
		}
	}
	else if (d[n][n] == 0) {
		for (int i = n; i < d.size(); i++) {
			genVect[i] = 1;
		}
		reverseRecGenVectorForZeroCase(d, p, n);
		return false;
	}
	else  {
		double first = sumOfDiagonalElements(d, n - 1);
		double second = -d[n][n];
		Matrix b = p.transpost();
		genVect[n] = sqrt(first / second);
		for (int i = n-2; i >= 0; i--) {
			double x = 0;
			for (int k = 1; k < n; k++) {
				x += b[i][k] * genVect[k];
			}
			genVect[i] = x;
		}
		//found solution
		return false;
	}
}





int main()
{
	vector<Matrix> startMatrix;
	vector<Matrix> A;
	vector<Matrix> B;
	vector<Matrix> C;

	for (int i = 0; i <= MAXN; i++) {
		int n = i, m = MAXN - i;
		generateMatrix(n, m, 0, &startMatrix);
	}
	/*
	Matrix firstMatrix = findMinMatrix(&startMatrix);
	std::cout << "Matrix with minimal determinant is: " << endl << firstMatrix << endl;
	
	generateVectorsFor(firstMatrix);
	for (Matrix vect : generatedVectors) {
		double result = matMult(vect, firstMatrix);
		if (result > 1) {
			A.push_back(firstMatrix);
		}
		else if (result == 1) {
			C.push_back(firstMatrix);
		}
		else if (result < 1) {
			B.push_back(firstMatrix);
		}
	}

	if (B.empty()) {
		cout << "First try" << endl;
		cout << "Min det for N = " << N << " is: " << firstMatrix.det() << endl;
		return 0;
	}


	A.clear();
	B.clear();
	C.clear();
	*/
	while (true) {
		for (Matrix matr : startMatrix) {
			//if (matr == firstMatrix) continue; // bad
			generateVectorsFor(matr);
			for (int i = 0; i < generatedVectors.size(); i++) {
				cout << "Generated vector i = " << endl << generatedVectors[i] << endl << endl;
			}
			for (Matrix vect : generatedVectors) {
				double result = matMult(vect, matr);
				if (result > 1) {
					A.push_back(matr);
				} else if (result == 1) {
					C.push_back(matr);
				} else if (result < 1) {
					B.push_back(matr);
				}
			}
		}

		if (B.empty()) {
			//если список матриц кончился, то все нашли
			break;
		}
		else {
			for (Matrix vect : generatedVectors) {
				for (Matrix a : A) {
					for (Matrix b : B) {
						Matrix tmp = calculateNewMatrix(a, b, vect);
						std::cout << endl << "New matrix = " << endl << tmp;
						C.push_back(tmp);
					}
				}
			}
			startMatrix.clear();
			startMatrix.assign(A.begin(), A.end()); //Возможно бага
			startMatrix.assign(C.begin(), C.end());
			A.clear();
			B.clear();
			C.clear();
		}
	}
	Matrix min = findMinMatrix(&startMatrix);
	std::cout << "Min det for N = " << N << " is: " << min.det() << endl;
	cout << "For matrix:" << endl << min << endl << endl;

	cout << "All other matrix determinants" << endl;
	for (Matrix m : startMatrix) {
		cout << m.det() << endl;
	}
	return 0;
}
