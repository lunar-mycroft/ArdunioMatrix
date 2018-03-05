#include"matrix.hpp"

bool canAdd(const Matrix & left, const Matrix & right) {
	return left.h_ == right.h_ && left.w_ == right.w_;
}

bool canMul(const Matrix & left, const Matrix & right) {
	return right.h == left.w;
}

Matrix::Matrix(float* in, unsigned char w, unsigned char h) : h_(h), w_(w) {
	for (unsigned char i = 0; i < w*h; ++i) {
		elements[i] = in[i];
	}
}

Matrix::Matrix(unsigned char w, unsigned char h) : h_(h), w_(w) {
	if (h*w == 0) {
		elements = nullptr;
	} else {
		elements = new float[h*w];
	}
};

Matrix::Matrix() : Matrix(0, 0) {};

Matrix::Matrix(const Matrix& other) {
	h_ = other.h_;
	w_ = other.w_;
	for (unsigned char i = 0; i < h_*w_; ++i) {
		elements[i] = other.elements[i];
	}
}

Matrix & Matrix::operator=(Matrix other) {
	swap(*this, other);
}

Matrix::~Matrix() {
	delete[] elements;
}

float & Matrix::operator[](loc coord) {
	return elements[coord.x, coord.y];
}

float Matrix::operator[](loc coord) const {
	return elements[coord.x, coord.y];
}

Matrix Matrix::operator+(const Matrix& other) const {
	Matrix result(w_, h_);
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_;++i) {
			result.elements[i] = other.elements[i]+elements[i];
		}
	} else {
		result.NaNify();
	}
}

Matrix & Matrix::operator+=(const Matrix& other) {
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_; ++i) {
			elements[i] = other.elements[i] + elements[i];
		}
	}
	else {
		NaNify();
	}
	return *this;
}

Matrix Matrix::operator-(const Matrix& other) const {
	Matrix result(w_, h_);
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_; ++i) {
			result.elements[i] = other.elements[i] - elements[i];
		}
	}
	else {
		result.NaNify();
	}
}

Matrix & Matrix::operator-=(const Matrix& other) {
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_; ++i) {
			elements[i] = other.elements[i] - elements[i];
		}
	}
	else {
		NaNify();
	}
	return *this;
}

Matrix Matrix::operator*(const float &coeff) const {
	Matrix result(w_, h_);
	for (unsigned char i = 0; i < h_*w_; ++i) {
		result.elements[i] = elements[i] * coeff;
	}
}

Matrix & Matrix::operator*=(const float& coeff) {
	for (unsigned char i = 0; i < h_*w_; ++i) {
		elements[i] *= coeff;
	}
}

Matrix Matrix::operator*(const Matrix & other) const {
	// from https://github.com/eecharlie/MatrixMath/blob/master/MatrixMath.cpp
	Matrix result(h_, other.w_);
	unsigned char i, j, k;
	for (i = 0; i < h_; i++) {
		for (j = 0; j < other.w_; j++)
		{
			result.elements[other.w_ * i + j] = 0;
			for (k = 0; k < w_; k++)
				result.elements[other.w_ * i + j] = result.elements[other.w_ * i + j] + elements[w_ * i + k] * other.elements[other.w_ * k + j];
		}
	}
}

Matrix & Matrix::operator*=(const Matrix & other) {
	(*this) = (*this)*other;
	return *this;
}

Matrix & Matrix::operator/=(const float& f) {
	(*this) = (*this) / f;
	return *this;
}