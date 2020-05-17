// ErmitConstant2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Matrix.h"
#include <set>
#include <algorithm>

using namespace std;

const int N = 3;
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
	double det = 10000;
	for (Matrix m : *matrixs) {
		if ((m.det() < det) && (m.det() > 0)) {
			cout << "m.det() = " << m.det() << endl;
			min = m;
			det = m.det();
		}
	}
	return min;
}

int findMinMatrixIndex(vector<Matrix>* matrixs) {
	Matrix min;
	int index = 0;
	for (int i = 0; i < matrixs->size(); i++) {
		Matrix m = matrixs->at(i);
		if (min.isNull()) {
			min = m;
		}
		else {
			if ((m.det() < min.det())) {
				min = m;
				index = i;
			}
		}
	}
	return index;
}

double matMult(Matrix vect, Matrix matr) {
	Matrix m1 = (vect.transpost() * matr) * vect;
	return m1.det();
}

Matrix calculateNewMatrix(Matrix f, Matrix h, Matrix vector) {

	double fx = matMult(vector, f);
	double hx = matMult(vector, h);
	cout << vector << endl;
	//fx and hx = 0
	if (fx == 0 && hx == 0) {
		return f; // bad
	}
	double c1 = 1 / (fx - hx);
	double c2 = 1 - hx;
	double c3 = fx - 1;
	Matrix tmp = c1 * ((c2 * f) + (c3 * h));
	return tmp;
}

bool isAllZerosVector(vector<int> v) {
	for (int i : v) {
		if (i != 0) return false;
	}
	return true;
}

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
			double coof = p[i][n - 1] * (-1);
			d = d.comb(i, n - 1, coof);
			p = p.comb(i, n - 1, coof);
		}
		d.updateDiagonalMatrix();
		//посчитать элемент вектора n-1
		double coof = d[n - 1][n - 1];
		double xi = 0;
		for (int i = n; i < d.size(); i++) {
			xi -= d[i][n - 1] / coof;
		}
		genVect[n - 1] = (int)xi;
		//рекурсивный вызов
		reverseRecGenVectorForZeroCase(d, p, n - 1);
	}
}

void higherZeroCaseRecFunc(Matrix d, Matrix p, int n, vector<int> vReverce, vector<Matrix>* vm) {
	if (n > 0) {
		//востановить
		//cout << "d.size - n-1 = " << d.size() - n - 1 << endl;
		//cout << "d" << endl << d << endl;
		//cout << "p" << endl << p << endl;
		//bad?
		for (int i = n + 1; i < d.size(); i++) {
			double coof = p[i][n] * (-1);
			d = d.comb(i, n, coof);
			p = p.comb(i, n, coof);
		}
		//

		d.updateDiagonalMatrix();
		//посчитать элемент вектора n-1

		//
		//cout << "d new" << endl << d << endl;
		//cout << "p new" << endl << p << endl;
		double a = 0; //bad

		for (int i = n + 1; i < d.size(); i++) {
			a += d[i][n] * v[i];
		}
		//
		double b = 1 - history;
		double c = p[n][n]; //bad?
		double xi = sqrt(abs(b / c)) - a;
		//cout << "history = " << history << endl;
		//cout << "xi = " << xi << endl;
		genVect[n] = abs((int)xi);
		double coof = d[n][n];
		history += coof * (pow(v[n] + a, 2));
		//рекурсивный вызов
		for (int k = (-1) * genVect[n]; k <= genVect[n]; k++) {
			vReverce.push_back(k);
			v[n] = k;
			//cout << "if (n > 0) k = " << k << ", v[" << n << "] = " << v[n] << endl;
			higherZeroCaseRecFunc(d, p, n - 1, vReverce, vm);
			vReverce.pop_back();
		}
	}
	else {
		for (int k = (-1) * genVect[n]; k <= genVect[n]; k++) {
			vReverce.push_back(k);
			v[n] = k;
			//cout << "else k = " << k << ", v[" << n << "] = " << v[n] << endl;
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
			//cout << "vm k = " << k << ", v[" << n << "] = " << v[n] << endl;
			higherZeroCaseRecFunc(d, p, n - 1, vReverce, &vm);
			vReverce.pop_back();

		}
		generatedVectors = vm;
	}
}

bool recGenerateVectors(Matrix d, Matrix p, int n) {
	if (d[n][n] > 0) {
		if (n < d.size() - 1) {
			//преобразования....
			for (int i = n + 1; i < d.size(); i++) {
				double coof = (d[i][n] / d[n][n]) * (-1);
				d = d.comb(i, n, coof);
				p = p.comb(i, n, coof);
			}
			//
			d.updateDiagonalMatrix();
			return recGenerateVectors(d, p, n + 1);
		}
		else {
			genVect[n] = abs(1 / (double)d[n][n]);
			v[n] = genVect[n];
			//cout << "1 / (double)d[n][n] = " << 1 / (double)d[n][n] << " v["<< n << "] = " << v[n] << endl;
			history = pow(genVect[n], 2) * d[n][n];
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
	else {
		double first = sumOfDiagonalElements(d, n - 1);
		double second = -d[n][n];
		Matrix b = p.transpost();
		genVect[n] = sqrt(first / second);
		for (int i = n - 2; i >= 0; i--) {
			double x = 0;
			for (int k = 1; k < n; k++) {
				x += b[k][i] * genVect[k];
			}
			genVect[i] = x;
		}
		//found solution
		return false;
	}
}


void cleanZerosVectors() {
	for (int i = 0; i < generatedVectors.size(); i++) {
		if (generatedVectors[i].isZeros()) {
			generatedVectors.erase(generatedVectors.begin() + i);
		}
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
	//sort matrix that first matrix will be with min det
	int minIndex = findMinMatrixIndex(&startMatrix);
	auto it = startMatrix.begin() + minIndex;
	rotate(startMatrix.begin(), startMatrix.begin() + minIndex, startMatrix.end());
	//startMatrix.swap(startMatrix[minIndex], startMatrix.begin()-1);

	for (int i = 0; i < startMatrix.size(); i++) {
		if (i > 15) {
			break;
		}
		Matrix matr = startMatrix[i];
		generateVectorsFor(matr);
		std::cout << "genereted new vetrors for i " << i << endl;
		std::cout << "Start matrix size = " << startMatrix.size() << endl;
		cleanZerosVectors();
		for (int i = 0; i < generatedVectors.size(); i++) {
			std::cout << "Generated vector i = " << endl << generatedVectors[i] << endl << endl;
		}
		/*
		if (generatedVectors.size() > 0) {
			std::cout << "B matrix = " << endl << matr << endl;
			//std::cout << "For vector = " << endl << vect << endl;
			B.push_back(matr);
		}
		*/

		//change this
		for (Matrix matr2 : startMatrix) {
			for (Matrix vect : generatedVectors) {
				double result = matMult(vect, matr2);
				if (result > 1) {
					A.push_back(matr2);
				}
				else if (result == 1) {
					C.push_back(matr2);
				}
				else if (result < 1) {
					std::cout << "B matrix = " << endl << matr << endl;
					std::cout << "For vector = " << endl << vect << endl;
					B.push_back(matr2);
				}
			}
		}

		//


		if (B.size() == 0) {
			//A.clear();
			//B.clear();
			//C.clear();
			//cout << "startMatrix size = " << endl;
			//если список матриц кончился, то все нашли
			//found = true;
			continue;
		}
		else {
			for (Matrix vect : generatedVectors) {
				for (Matrix a : A) {
					for (Matrix b : B) {
						Matrix tmp = calculateNewMatrix(a, b, vect);
						std::cout << endl << "New matrix = " << endl << tmp << endl;
						std::cout << endl << "For matrix A = " << endl << a << endl;
						std::cout << endl << "And B = " << endl << b << endl;
						C.push_back(tmp);
					}
				}
			}
			startMatrix.insert(startMatrix.end(), A.begin(), A.end());
			startMatrix.insert(startMatrix.end(), C.begin(), C.end());
			A.clear();
			B.clear();
			C.clear();
		}
	}
	Matrix min = findMinMatrix(&startMatrix);
	std::cout << "Min det for N = " << N << " is: " << min.det() << endl;
	std::cout << "For matrix:" << endl << min << endl << "with size " << min.size() << endl << endl;

	std::cout << "All other matrix determinants" << endl;
	for (Matrix m : startMatrix) {
		std::cout << m.det() << endl;
	}
	return 0;
}
