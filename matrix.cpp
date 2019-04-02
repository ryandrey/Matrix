#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

template <typename T>
class Matrix;

template <typename T>
class MatrixIterator {
    private:
        Matrix<T>& matrix;
        size_t row, column;

    public:
        MatrixIterator(Matrix<T>& m, size_t i, size_t j)
            : matrix(m)
            , row(i)
            , column(j) {
        }

        bool operator == (MatrixIterator other) const {
            return row == other.row
                   && column == other.column;
        }

        bool operator != (MatrixIterator other) const {
            return !(*this == other);
        }

        T& operator * () {
            return matrix[row][column];
        }

        MatrixIterator& operator++ () {
            ++column;
            if (column == matrix.size().second) {
                ++row;
                column = 0;
            }
            return *this;
        }
};

template <typename T>
class ConstMatrixIterator {
private:
    const Matrix<T>& matrix;
    size_t row, column;

public:
    ConstMatrixIterator(const Matrix<T>& m, size_t i, size_t j)
            : matrix(m)
            , row(i)
            , column(j) {
    }

    bool operator == (ConstMatrixIterator other) const {
        return row == other.row
               && column == other.column;
    }

    bool operator != (ConstMatrixIterator other) const {
        return !(*this == other);
    }

    const T& operator * () {
        return matrix[row][column];
    }

    ConstMatrixIterator& operator++ () {
        ++column;
        if (column == matrix.size().second) {
            ++row;
            column = 0;
        }
        return *this;
    }
};

template <typename T>
class Matrix {
    private:
        vector<vector<T>> data;

    public:
        Matrix(size_t rows, size_t columns) {
            data = vector<vector<T>>(rows, vector<T>(columns, 0));
        }

        MatrixIterator<T> begin() {
            return {*this, 0, 0};
        }

        MatrixIterator<T> end() {
            return {*this, this->size().first, 0};
        }

        ConstMatrixIterator<T> begin() const {
            return {*this, 0, 0};
        }

        ConstMatrixIterator<T> end() const {
            return {*this, this->size().first, 0};
        }

        Matrix(const vector<vector<T>>& d)
            : data(d) {
        }

        const vector<T>& operator[] (size_t row) const {
            return data[row];
        }

        vector<T>& operator[] (size_t row) {
            return data[row];
        }

        const pair<size_t, size_t> size() const {
            return {data.size(), data[0].size()};
        }

        Matrix<T> operator + (const Matrix<T>& other) const {
            vector<vector<T>> v(other.size().first, vector<T>(other.size().second, 0));
            Matrix<T> result(v);
            for (size_t i = 0; i != this->size().first; ++i) {
                for (size_t j = 0; j != this->size().second; ++j) {
                    result[i][j] = data[i][j] + other[i][j];
                }
            }
            return result;
        }

        Matrix<T>& operator += (const Matrix<T>& other) {
            for (size_t i = 0; i != this->size().first; ++i) {
                for (size_t j = 0; j != this->size().second; ++j) {
                    (*this)[i][j] += other[i][j];
                }
            }
            return *this;
        }

        Matrix<T>& transpose() {
            vector<vector<T>> v(this->size().second, vector<T>(this->size().first, 0));
            for (size_t i = 0; i != this->size().first; ++i) {
                for (size_t j = 0; j != this->size().second; ++j) {
                    v[j][i] = (*this)[i][j];
                }
            }
            *this = Matrix<T>(v);
            return *this;
        }

        Matrix<T> transposed() const {
            vector<vector<T>> v(this->size().second, vector<T>(this->size().first, 0));
            for (size_t i = 0; i != this->size().first; ++i) {
                for (size_t j = 0; j != this->size().second; ++j) {
                    v[j][i] = (*this)[i][j];
                }
            }
            Matrix<T> result(v);
            return result;
        }

        Matrix<T>& operator *= (const Matrix<T>& other) {
            assert(other.size().first == this->size().second);
            vector<vector<T>> v(this->size().first, vector<T>(other.size().second, 0));
            Matrix<T> result(v);
            for (size_t i = 0; i != this->size().first; ++i) {
                for (size_t j = 0; j != other.size().second; ++j) {
                    for (size_t t = 0; t != other.size().first; ++t) {
                        result[i][j] += (*this)[i][t] * other[t][j];
                    }
                }
            }
            *this = result;
            return *this;
        }

    Matrix<T> operator * (const Matrix<T>& other) const {
        assert(other.size().first == this->size().second);
        vector<vector<T>> v(this->size().first, vector<T>(other.size().second, 0));
        Matrix<T> result(v);
        for (size_t i = 0; i != this->size().first; ++i) {
            for (size_t j = 0; j != other.size().second; ++j) {
                for (size_t t = 0; t != other.size().first; ++t) {
                    result[i][j] += (*this)[i][t] * other[t][j];
                }
            }
        }
        return result;
    }
};

template <typename T, typename Scalar>
Matrix<T>& operator *= (Matrix<T>& m, Scalar scalar) {
    for (size_t i = 0; i != m.size().first; ++i) {
        for (size_t j = 0; j != m.size().second; ++j) {
            m[i][j] *= static_cast<T>(scalar);
        }
    }
    return m;
}

template <typename T, typename Scalar>
Matrix<T> operator * (const Matrix<T>& m, Scalar scalar) {
    vector<vector<T>> v(m.size().first, vector<T>(m.size().second, 0));
    Matrix<T> result(v);
    for (size_t i = 0; i != m.size().first; ++i) {
        for (size_t j = 0; j != m.size().second; ++j) {
            result[i][j] = m[i][j] * static_cast<T>(scalar);
        }
    }
    return result;
}

template <typename T>
ostream& operator << (ostream& out, const Matrix<T> m) {
    for (size_t i = 0; i != m.size().first; ++i) {
        for (size_t j = 0; j != m.size().second; ++j) {
            if (j > 0) {
                out << '\t';
            }
            out << m[i][j];
        }
        if (i < m.size().first - 1) {
            out << '\n';
        }
    }
    return out;
}
