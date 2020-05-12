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
	Matrix(int rows, int cols, const double& init = double());

	Matrix(int dim, double* aij);
	Matrix(vector<int> v);

	Matrix(int rows, int cols, double** matrix);
	Matrix(const Matrix& m);
	bool isNull();
	bool isSquare();
	bool isDiagonal();
	bool isIdentity();
	bool isStep();
	bool isRref();
	bool isInvertible();

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

/* CONSTRUCTORS */

Matrix::Matrix(int dim, const double& init) {
	if (dim < 1) {
		cerr << "[Matrix(int, const double&)] invalid argument exception: dim must be positive" << endl;
		throw "[Matrix(int, const double&)] invalid argument exception: dim must be positive";
	}
	this->rows = dim;
	this->cols = dim;
	this->matrix = new double * [dim];
	for (int i = 0; i < dim; i++) {
		this->matrix[i] = new double[dim];
		for (int j = 0; j < dim; j++) {
			this->matrix[i][j] = init;
		}
	}
}


Matrix::Matrix(int dim, double* aij) {
	if (dim < 1) {
		cerr << "[Matrix(int, const double&)] invalid argument exception: dim must be positive" << endl;
		throw "[Matrix(int, const double&)] invalid argument exception: dim must be positive";
	}
	int counter = 0;
	this->rows = dim;
	this->cols = dim;
	this->matrix = new double* [dim];
	for (int i = 0; i < dim; i++) {
		this->matrix[i] = new double[dim];
	}

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j <= i; j++) {
			if (i == j) this->matrix[i][j] = 1;
			else { 
				this->matrix[i][j] = aij[counter];
				this->matrix[j][i] = aij[counter];
				counter++;
			}
		}
	}
}


Matrix::Matrix(vector<int> v){
	int rows = v.size();
	if (rows < 1) {
		cerr << "[Matrix(int, const double&)] invalid argument exception: rows must be positive" << endl;
		throw "[Matrix(int, const double&)] invalid argument exception: rows must be positive";
	}
	this->rows = rows;
	this->cols = 1;
	this->matrix = new double* [rows];
	for (int i = 0; i < rows; i++) {
		this->matrix[i] = new double[cols];
		for (int j = 0; j < cols; j++) {
			this->matrix[i][j] = v[i];
		}
	}

}

Matrix::Matrix(int rows, int cols, const double& init) {
	if (rows < 1) {
		cerr << "[Matrix(int, int, const double&)] invalid argument exception: rows must be positive" << endl;
		throw "[Matrix(int, int, const double&)] invalid argument exception: rows must be positive";
	}
	if (cols < 1) {
		cerr << "[Matrix(int, int, const double&)] invalid argument exception: cols must be positive" << endl;
		throw "[Matrix(int, int, const double&)] invalid argument exception: cols must be positive";
	}
	this->rows = rows;
	this->cols = cols;
	this->matrix = new double * [rows];
	for (int i = 0; i < rows; i++) {
		this->matrix[i] = new double[cols];
		for (int j = 0; j < cols; j++) {
			this->matrix[i][j] = init;
		}
	}
}

Matrix::Matrix(int rows, int cols, double** matrix) {
	if (rows < 1) {
		cerr << "[Matrix(int, int, const double&)] invalid argument exception: rows must be positive" << endl;
		throw "[Matrix(int, int, const double&)] invalid argument exception: rows must be positive";
	}
	if (cols < 1) {
		cerr << "[Matrix(int, int, const double&)] invalid argument exception: cols must be positive" << endl;
		throw "[Matrix(int, int, const double&)] invalid argument exception: cols must be positive";
	}
	this->rows = rows;
	this->cols = cols;
	this->matrix = new double * [rows];
	for (int i = 0; i < this->rows; i++) {
		this->matrix[i] = new double[cols];
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j] = matrix[i][j];
		}
	}
}

Matrix::Matrix(const Matrix& m) {
	const int rows = m.rows;
	const int cols = m.cols;
	this->rows = rows;
	this->cols = cols;
	this->matrix = new double * [rows];
	for (int i = 0; i < this->rows; i++) {
		this->matrix[i] = new double[cols];
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j] = m.matrix[i][j];
		}
	}
}

/* PRIVATE METHODS */

void Matrix::_mStep(int row, int col) {
	if (row < this->rows - 1 && col < this->cols) {
		int ref = -1;
		for (int i = row; i < this->rows && ref == -1; i++) {
			if (this->matrix[i][col] != 0) {
				ref = i;
			}
		}
		if (ref != -1) {
			if (row != ref) {
				this->mSwap(row, ref);
			}
			if (ref != this->rows - 1) {
				for (int i = row + 1; i < this->rows; i++) {
					if (this->matrix[i][col] != 0) {
						this->mComb(i, row, -this->matrix[i][col] / this->matrix[row][col]);
					}
				}
			}
			this->_mStep(row + 1, col + 1);
		}
		else {
			this->_mStep(row, col + 1);
		}
	}
}

void Matrix::_mRref(int r, int c) {
	this->mStep();
	const int row = this->rows - r - 1;
	const int col = this->cols - c - 1;
	if (col >= 0 && row >= 0) {
		int ref = -1;
		for (int i = 0; i <= col && ref == -1; i++) {
			if (this->matrix[row][i] != 0) {
				ref = i;
			}
		}
		if (ref == -1) {
			this->_mRref(r + 1, c);
		}
		else {
			if (this->matrix[row][ref] != 1) {
				this->mMult(row, 1 / this->matrix[row][ref]);
			}
			for (int i = 0; i < row; i++) {
				if (this->matrix[i][ref] != 0) {
					this->mComb(i, row, -this->matrix[i][ref]);
				}
			}
			this->_mRref(r + 1, c + col - ref + 1);
		}
	}
}

void Matrix::_iStep(Matrix& inverse, int row, int col) {
	if (row < this->rows - 1 && col < this->cols) {
		int ref = -1;
		for (int i = row; i < this->rows && ref == -1; i++) {
			if (this->matrix[i][col] != 0) {
				ref = i;
			}
		}
		if (ref != -1) {
			this->mSwap(row, ref);
			inverse.mSwap(row, ref);
			if (ref != this->rows - 1) {
				for (int i = row + 1; i < this->rows; i++) {
					if (this->matrix[i][col] != 0) {
						const double k = -this->matrix[i][col] / this->matrix[row][col];
						this->mComb(i, row, k);
						inverse.mComb(i, row, k);
					}
				}
			}
			this->_iStep(inverse, row + 1, col + 1);
		}
		else {
			this->_iStep(inverse, row, col + 1);
		}
	}
}

void Matrix::_iRref(Matrix& inverse, int r, int c) {
	const int row = this->rows - r - 1;
	const int col = this->cols - c - 1;
	if (col >= 0 && row >= 0) {
		int ref = -1;
		for (int i = 0; i <= col && ref == -1; i++) {
			if (this->matrix[row][i] != 0) {
				ref = i;
			}
		}
		if (ref == -1) {
			this->_iRref(inverse, r + 1, c);
		}
		else {
			if (this->matrix[row][ref] != 1) {
				const double k = 1 / this->matrix[row][ref];
				this->mMult(row, k);
				inverse.mMult(row, k);
			}
			for (int i = 0; i < row; i++) {
				if (this->matrix[i][ref] != 0) {
					const double k = -this->matrix[i][ref];
					this->mComb(i, row, k);
					inverse.mComb(i, row, k);
				}
			}
			this->_iRref(inverse, r + 1, c + col - ref + 1);
		}
	}
}

double Matrix::_det(Matrix& m) {
	double det = 0;
	if (m.rows == 1 && m.cols == 1) {
		det = m.matrix[0][0];
	}
	else if (m.rows == 2 && m.cols == 2) {
		det = m.matrix[0][0] * m.matrix[1][1] - m.matrix[0][1] * m.matrix[1][0];
	}
	else {
		bool isRow = rand() % 2;
		if (isRow) {
			int row = rand() % m.rows;
			for (int i = 0; i < m.cols; i++) {
				if (m.matrix[row][i] != 0) {
					Matrix s = m.sub(row, i);
					int n = (i + row) % 2 == 0 ? 1 : -1;
					det += n * m.matrix[row][i] * m._det(s);
				}
			}
		}
		else {
			int col = rand() % m.cols;
			for (int i = 0; i < m.rows; i++) {
				if (m.matrix[i][col] != 0) {
					Matrix s = m.sub(i, col);
					int n = (i + col) % 2 == 0 ? 1 : -1;
					det += n * m.matrix[i][col] * m._det(s);
				}
			}
		}
	}
	return det;
}

/* PUBLIC METHODS */

int Matrix::size() {
	return cols;
}

void Matrix::print() {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			cout << this->matrix[i][j] << " ";
		}
		cout << endl;
	}
}

bool Matrix::isNull() {
	bool res = true;
	for (int i = 0; i < this->rows && res; i++) {
		for (int j = 0; j < this->cols && res; j++) {
			if (this->matrix[i][j] != 0) {
				res = false;
			}
		}
	}
	return res;
}

bool Matrix::isSquare() {
	return (this->rows == this->cols);
}

bool Matrix::isDiagonal() {
	bool res = true;
	if (!this->isSquare()) {
		res = false;
	}
	for (int i = 0; i < this->rows && res; i++) {
		for (int j = 0; j < this->cols && res; j++) {
			if (i != j && this->matrix[i][j] != 0 || i == j && this->matrix[i][j] == 0) {
				res = false;
			}
		}
	}
	return res;
}

bool Matrix::isIdentity() {
	bool res = true;
	if (!this->isSquare()) {
		res = false;
	}
	for (int i = 0; i < this->rows && res; i++) {
		for (int j = 0; j < this->cols && res; j++) {
			if ((i != j && this->matrix[i][j] != 0) || (i == j && this->matrix[i][j] != 1)) {
				res = false;
			}
		}
	}
	return res;
}

bool Matrix::isStep() {
	bool res = true;
	for (int i = 0, col = -1; i < this->rows && res; i++) {
		bool stop = false;
		for (int j = 0; j < this->cols && !stop; j++) {
			if (this->matrix[i][j] != 0) {
				if (j <= col) {
					res = false;
				}
				col = j;
				stop = true;
			}
			else if (j == this->cols - 1) {
				col = this->cols;
			}
		}
	}
	return res;
}

bool Matrix::isRref() {
	bool res = true;
	for (int i = 0, col = -1; i < this->rows && res; i++) {
		bool pivot = false;
		for (int j = 0; j < this->cols && res; j++) {
			if (this->matrix[i][j] != 0) {
				if (j < col || this->matrix[i][j] != 1 || pivot) {
					res = false;
				}
				else {
					col = j;
					pivot = true;
				}

			}
			else if (j == this->cols - 1 && !pivot) {
				col = this->cols;
			}
		}
	}
	return res;
}

bool Matrix::isInvertible() {
	return (this->isSquare() && this->rang() == this->rows);
}

void Matrix::mSwap(int x, int y) {
	if (x >= this->rows || y >= this->rows) {
		cerr << "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows" << endl;
		throw "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows";
	}
	else {
		double* temp = this->matrix[x];
		this->matrix[x] = this->matrix[y];
		this->matrix[y] = temp;
	}
}

void Matrix::mMult(int x, double k) {
	if (x >= this->rows) {
		cerr << "[Matrix::mSwap(int, int)] invalid agrument exceptions: x must be lower than rows" << endl;
		throw "[Matrix::mSwap(int, int)] invalid agrument exceptions: x must be lower than rows";
	}
	else {
		for (int i = 0; i < this->cols; i++) {
			this->matrix[x][i] *= k;
		}
	}
}

void Matrix::mComb(int x, int y, double k) {
	if (x >= this->rows || y >= this->rows) {
		cerr << "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows" << endl;
		throw "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows";
	}
	else {
		for (int i = 0; i < this->cols; i++) {
			this->matrix[x][i] += k * this->matrix[y][i];
		}
	}
}

Matrix Matrix::swap(int x, int y) {
	if (x >= this->rows || y >= this->rows) {
		cerr << "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows" << endl;
		throw "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows";
	}
	else {
		Matrix m = *this;
		m.mSwap(x, y);
		return m;
	}
}

Matrix Matrix::mult(int x, double k) {
	if (x >= this->rows) {
		cerr << "[Matrix::mSwap(int, int)] invalid agrument exceptions: x must be lower than rows" << endl;
		throw "[Matrix::mSwap(int, int)] invalid agrument exceptions: x must be lower than rows";
	}
	else {
		Matrix m = *this;
		m.mMult(x, k);
		return m;
	}
}

Matrix Matrix::comb(int x, int y, double k) {
	if (x >= this->rows || y >= this->rows) {
		cerr << "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows" << endl;
		throw "[Matrix::mSwap(int, int)] invalid agrument exceptions: x and y must be lower than rows";
	}
	else {
		Matrix m = *this;
		m.mComb(x, y, k);
		return m;
	}
}

Matrix Matrix::transpost() {
	double** m = new double * [this->cols];
	for (int i = 0; i < this->cols; i++) {
		m[i] = new double[this->rows];
		for (int j = 0; j < this->rows; j++) {
			m[i][j] = this->matrix[j][i];
		}
	}
	return Matrix(this->cols, this->rows, m);
}

void Matrix::updateDiagonalMatrix() {
	for (int i = 0; i < this->cols; i++) {
		for (int j = 0; j <= i; j++) {
			if (i != j) {
				this->matrix[j][i] = this->matrix[i][j];
			}
		}
	}
}


void Matrix::mStep() {
	this->_mStep();
}

Matrix Matrix::step() {
	Matrix m = *this;
	m.mStep();
	return m;
}

void Matrix::mRref() {
	this->_mRref();
}

Matrix Matrix::rref() {
	Matrix m = *this;
	m.mRref();
	return m;
}

Matrix Matrix::mInvert() {
	if (!this->isInvertible()) {
		cerr << "[Matrix mInvert] invalid matrix exception: matrix must be invertible" << endl;
		throw "[Matrix mInvert] invalid matrix exception: matrix must be invertible";
	}
	Matrix inverse = identity(this->rows);
	this->_iStep(inverse);
	this->_iRref(inverse);
	return inverse;
}

Matrix Matrix::invert() {
	Matrix m = *this;
	if (!m.isInvertible()) {
		cerr << "[Matrix invert] invalid matrix exception: matrix must be invertible" << endl;
		throw "[Matrix invert] invalid matrix exception: matrix must be invertible";
	}
	Matrix inverse = identity(m.rows);
	m._iStep(inverse);
	m._iRref(inverse);
	return inverse;
}

int Matrix::rang() {
	Matrix m = *this;
	int rang = 0;
	m.mStep();
	for (int i = 0; i < m.rows; i++) {
		bool stop = false;
		for (int j = 0; j < m.cols && !stop; j++) {
			if (m.matrix[i][j] != 0) {
				rang++;
				stop = true;;
			}
		}
	}
	return rang;
}

Matrix Matrix::sub(int row, int col) {
	if (this->rows == 1 || this->cols == 1) {
		cerr << "[Matrix sub] invalid matrix exception: matrix must have more than one rows and cols" << endl;
		throw "[Matrix sub] invalid matrix exception: matrix must have more than one rows and cols";
	}
	if (row < 0 || col < 0 || row >= this->rows || col >= this->cols) {
		cerr << "[Matrix sub] invalid argument exception: invalid row or col" << endl;
		throw "[Matrix sub] invalid argument exception: invalid row or col";
	}
	double** m = new double * [this->rows - 1];
	for (int i = 0, k = 0; i < this->rows; i++) {
		if (i != row) {
			m[k] = new double[this->cols - 1];
			for (int j = 0, l = 0; j < this->cols; j++) {
				if (j != col) {
					m[k][l] = this->matrix[i][j];
					l++;
				}
			}
			k++;
		}
	}
	return Matrix(this->rows - 1, this->cols - 1, m);
}

double Matrix::det() {
	if (!this->isSquare()) {
		cerr << "[Matrix det] invalid matrix exception: matrix must be square" << endl;
		throw "[Matrix det] invalid matrix exception: matrix must be square";
	}
	srand(time(NULL));
	return this->_det(*this);
}

Matrix& Matrix::operator=(const Matrix& m) {
	const int rows = m.rows;
	const int cols = m.cols;
	this->rows = rows;
	this->cols = cols;
	this->matrix = new double * [rows];
	for (int i = 0; i < this->rows; i++) {
		this->matrix[i] = new double[cols];
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j] = m.matrix[i][j];
		}
	}
	return *this;
}

Matrix Matrix::operator-() {
	double** m = new double * [this->rows];
	for (int i = 0; i < this->rows; i++) {
		m[i] = new double[this->cols];
		for (int j = 0; j < this->cols; j++) {
			m[i][j] = -this->matrix[i][j];
		}
	}
	return Matrix(this->rows, this->cols, m);
}

Matrix& Matrix::operator++() {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j]++;
		}
	}
	return *this;
}

Matrix& Matrix::operator--() {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j]--;
		}
	}
	return *this;
}

Matrix Matrix::operator++(int n) {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j]++;
		}
	}
	return *this;
}

Matrix Matrix::operator--(int n) {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j]--;
		}
	}
	return *this;
}

Matrix Matrix::operator+(const Matrix& matrix) {
	if (this->rows != matrix.rows || this->cols != matrix.cols) {
		cerr << "[Matrix operator+] invalid argument exception: rows and cols cannot be different" << endl;
		throw "[Matrix operator+] invalid argument exception: rows and cols cannot be different";
	}
	else {
		double** m = new double * [this->rows];
		for (int i = 0; i < this->rows; i++) {
			m[i] = new double[this->cols];
			for (int j = 0; j < this->cols; j++) {
				m[i][j] = this->matrix[i][j] + matrix.matrix[i][j];
			}
		}
		Matrix res = Matrix(this->rows, this->cols, m);
		return res;
	}
}

Matrix Matrix::operator-(const Matrix& matrix) {
	if (this->rows != matrix.rows || this->cols != matrix.cols) {
		cerr << "[Matrix operator-] invalid argument exception: rows and cols cannot be different" << endl;
		throw "[Matrix operator-] invalid argument exception: rows and cols cannot be different";
	}
	else {
		double** m = new double * [this->rows];
		for (int i = 0; i < this->rows; i++) {
			m[i] = new double[this->cols];
			for (int j = 0; j < this->cols; j++) {
				m[i][j] = this->matrix[i][j] - matrix.matrix[i][j];
			}
		}
		return Matrix(this->rows, this->cols, m);
	}
}

Matrix Matrix::operator*(const Matrix& matrix) {
	if (this->cols != matrix.rows) {
		cerr << "[Matrix operator*(const Matrix&)] invalid argument exception: not conformable matrixes" << endl;
		throw "[Matrix operator*(const Matrix&)] invalid argument exception: not conformable matrixes";
	}
	else {
		double** m = new double * [this->rows];
		for (int i = 0; i < this->rows; i++) {
			m[i] = new double[matrix.cols];
			for (int j = 0; j < matrix.cols; j++) {
				m[i][j] = 0;
				for (int k = 0; k < this->cols; k++) {
					m[i][j] += this->matrix[i][k] * matrix.matrix[k][j];
				}
			}
		}
		return Matrix(this->rows, matrix.cols, m);
	}
}

Matrix Matrix::operator*(double k) {
	double** m = new double * [this->rows];
	for (int i = 0; i < this->rows; i++) {
		m[i] = new double[this->cols];
		for (int j = 0; j < this->cols; j++) {
			m[i][j] = this->matrix[i][j] * k;
		}
	}
	return Matrix(this->rows, this->cols, m);
}

Matrix Matrix::operator/(double k) {
	double** m = new double * [this->rows];
	for (int i = 0; i < this->rows; i++) {
		m[i] = new double[this->cols];
		for (int j = 0; j < this->cols; j++) {
			m[i][j] = this->matrix[i][j] / k;
		}
	}
	return Matrix(this->rows, this->cols, m);
}

Matrix Matrix::operator|(const Matrix& matrix) {
	if (this->rows != matrix.rows) {
		cerr << "[Matrix::operator|] invalid agrument exceptions: matrix must have same rows number" << endl;
		throw "[Matrix::operator|] invalid agrument exceptions: matrix must have same rows number";
	}
	const int rows = this->rows;
	const int cols = this->cols + matrix.cols;
	double** m = new double * [rows];
	for (int i = 0; i < rows; i++) {
		m[i] = new double[cols];
		for (int j = 0; j < cols; j++) {
			if (this->cols - j > 0) {
				m[i][j] = this->matrix[i][j];
			}
			else {
				m[i][j] = matrix.matrix[i][j - this->cols];
			}
		}
	}
	return Matrix(rows, cols, m);
}

Matrix Matrix::operator+=(const Matrix& matrix) {
	if (this->rows != matrix.rows || this->cols != matrix.cols) {
		cerr << "[Matrix operator+=] invalid argument exception: rows and cols cannot be different" << endl;
		throw "[Matrix operator+=] invalid argument exception: rows and cols cannot be different";
	}
	else {
		for (int i = 0; i < this->rows; i++) {
			for (int j = 0; j < this->cols; j++) {
				this->matrix[i][j] += matrix.matrix[i][j];
			}
		}
		return *this;
	}
}

Matrix Matrix::operator-=(const Matrix& matrix) {
	if (this->rows != matrix.rows || this->cols != matrix.cols) {
		cerr << "[Matrix operator-=] invalid argument exception: rows and cols cannot be different" << endl;
		throw "[Matrix operator-=] invalid argument exception: rows and cols cannot be different";
	}
	else {
		for (int i = 0; i < this->rows; i++) {
			for (int j = 0; j < this->cols; j++) {
				this->matrix[i][j] -= matrix.matrix[i][j];
			}
		}
		return *this;
	}
}

Matrix Matrix::operator*=(double k) {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j] *= k;
		}
	}
	return *this;
}

Matrix Matrix::operator/=(double k) {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			this->matrix[i][j] /= k;
		}
	}
	return *this;
}

bool Matrix::operator==(const Matrix& matrix) {
	bool res = true;
	if (this->rows != matrix.rows || this->cols != matrix.cols) {
		res = false;
	}
	int i, j;
	for (i = 0; i < this->rows && res; i++) {
		for (j = 0; j < this->cols && res; j++) {
			if (this->matrix[i][j] != matrix.matrix[i][j]) {
				res = false;
			}
		}
	}
	return res;
}

bool Matrix::operator!=(const Matrix& matrix) {
	bool res = false;
	if (this->rows != matrix.rows || this->cols != matrix.cols) {
		res = true;
	}
	int i, j;
	for (i = 0; i < this->rows && !res; i++) {
		for (j = 0; j < this->cols && !res; j++) {
			if (this->matrix[i][j] != matrix.matrix[i][j]) {
				res = true;
			}
		}
	}
	return res;
}

//This is fix for error with set and loops
bool operator<(const Matrix& m, const Matrix& m2) {
	int cm = 0;
	int cm2 = 0;
	if (m.cols < m2.cols && m.rows < m2.rows) return true;
	else if (m.cols > m2.cols && m.rows > m2.rows) return false;
	else if (m.cols == m2.cols && m.rows == m2.rows) {
		for (int i = 0; i < m.rows; i++) {
			for (int j = 0; j < m.cols; j++) {
				if (m.matrix[i][j] < m2.matrix[i][j]) cm++;
				else if (m.matrix[i][j] > m2.matrix[i][j]) cm2++;
			}
		}
	}
	return cm < cm2;
}

double* Matrix::operator[] (int i) {
	if (i < 0 || i >= this->rows) {
		cerr << "[Matrix operator[]] invalid argument exception: invalid index" << endl;
		throw "[Matrix operator[]] invalid argument exception: invalid index";
	}
	return this->matrix[i];
}

ostream& operator<< (ostream& out, const Matrix& matrix) {
	for (int i = 0; i < matrix.rows; i++) {
		for (int j = 0; j < matrix.cols; j++) {
			out << matrix.matrix[i][j] << "\t";
		}
		out << endl;
	}
	return out;
}

istream& operator>> (istream& in, Matrix& matrix) {
	for (int i = 0; i < matrix.rows; i++) {
		for (int j = 0; j < matrix.cols; j++) {
			in >> matrix.matrix[i][j];
		}
	}
	return in;
}

Matrix operator*(double k, const Matrix& m) {
	double** matrix = new double * [m.rows];
	for (int i = 0; i < m.rows; i++) {
		matrix[i] = new double[m.cols];
		for (int j = 0; j < m.cols; j++) {
			matrix[i][j] = m.matrix[i][j] * k;
		}
	}
	return Matrix(m.rows, m.cols, matrix);
}

Matrix diagonal(int n, double k) {
	double** matrix = new double * [n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (int j = 0; j < n; j++) {
			matrix[i][j] = i == j ? k : 0;
		}
	}
	return Matrix(n, n, matrix);
}

Matrix identity(int n) {
	double** matrix = new double * [n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (int j = 0; j < n; j++) {
			matrix[i][j] = i == j ? 1 : 0;
		}
	}
	return Matrix(n, n, matrix);
}
