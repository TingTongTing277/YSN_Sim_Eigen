/*
 * Matrix.h
 *
 *  Created on: 2023��4��9��
 *      Author: fuzic
 */


#ifndef MATH_MATRIX_H_
#define MATH_MATRIX_H_

//#include "Math_Header.h"
#include <stdexcept>  // ����쳣����
#include <iostream>
#include <iomanip>
#include <string>

namespace Sim_Eigen {
    template<typename Scalar>
    class Block;

    template<typename Scalar>
    class Vector3;

    template<typename Scalar>
    class Matrix
    {
    private:
        int row;
        int column;
        Scalar* matrix;

        // ���ų�ʼ����ʵ��
        template<typename, typename> friend class CommaInitializer;

        template<typename Derived, typename Scalar2>
        class CommaInitializer {
            Derived& matrix_;
            int index_;
            int expectedSize_;
        public:
            CommaInitializer(Derived& mat, Scalar2 first_val) : matrix_(mat), index_(0) {
                expectedSize_ = mat.rows() * mat.cols();
                *this, first_val;  // �����һ��ֵ
            }

            CommaInitializer& operator,(Scalar2 value) {
                if (index_ >= expectedSize_) {
                    throw std::out_of_range("Too many elements during initialization");
                }
                matrix_.data()[index_++] = static_cast<Scalar2>(value);
                return *this;
            }

            ~CommaInitializer() noexcept(false) {
                if (index_ != expectedSize_) {
                    throw std::invalid_argument(
                        "Incorrect number of elements: " +
                        std::to_string(index_) + " given, " +
                        std::to_string(expectedSize_) + " expected");
                }
            }
        };
    public:
        Matrix(int x, int y) :row(x), column(y)//������
        {
            if (x > 0 && y > 0) {
                matrix = new Scalar[row * column];
                setZero();
            }
            else {
                matrix = nullptr;
            }
        }
        Matrix(Matrix& other) :row(other.row), column(other.column)//�������캯��
        {
            matrix = new Scalar[row * column];
            std::copy(other.matrix, other.matrix + row * column, matrix);
        }

        // �ƶ����캯��
        Matrix(Matrix&& other) noexcept : row(other.row), column(other.column), matrix(other.matrix){//, mp(other.mp) {
            other.matrix = nullptr;
            other.row = 0;
            other.column = 0;
        }

        ~Matrix()
        {
            delete[] matrix;
        }

    public:
        // Eigen ���ӿ�
        int rows() const { return row; }
        int cols() const { return column; }
        Scalar* data() { return matrix; }
        const Scalar* data() const { return matrix; }

        // �����Ԫ����
        friend std::ostream& operator<<(std::ostream& os, const Matrix<Scalar>& mat) {
            // ����ԭ�и�ʽ״̬
            std::ios_base::fmtflags origFlags = os.flags();
            std::streamsize origPrecision = os.precision();
            std::streamsize origWidth = os.width();

            // ���������ʽ��ʾ�����̶�С���㣬4λ���ȣ�
            os << std::fixed << std::setprecision(4);

            for (int i = 0; i < mat.rows(); ++i) {
                for (int j = 0; j < mat.cols(); ++j) {
                    os << mat(i, j);
                    if (j != mat.cols() - 1) os << "\t";  // Ԫ�ؼ��ÿո�ָ�
                }
                if (i != mat.rows() - 1) os << std::endl;  // ��ĩ����
            }
            os << std::endl; // ��ĩ����

            // �ָ�ԭ�и�ʽ
            os.flags(origFlags);
            os.precision(origPrecision);
            os.width(origWidth);
            return os;
        }

        // ���������
        // ��ʽ��ʼ�������
        CommaInitializer<Matrix<Scalar>, Scalar> operator<<(Scalar value) {
            if (!matrix) throw std::runtime_error("Matrix not initialized");
            return CommaInitializer<Matrix<Scalar>, Scalar>(*this, value);
        }

        //��������
        inline Scalar& operator()(int rows, int cols) {
            if (rows >= row || cols >= column || rows < 0 || cols < 0) {
                throw std::invalid_argument("Matrix element access is illegal");
            }
            return matrix[rows*column+cols];
        }

        inline const Scalar& operator()(int rows, int cols) const {
            if (rows >= row || cols >= column || rows < 0 || cols < 0) {
                throw std::invalid_argument("Matrix element access is illegal");
            }
            return matrix[rows*column+cols];
        }

        Matrix operator + (const Matrix& m) const
        {
            if (row != m.row || column != m.column) {
                throw std::invalid_argument("Matrix addition: Matrix dimensions mismatch");
            }
            Matrix r(row, column);
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    r.matrix[i*r.column + j] = m.matrix[i*m.column + j] + matrix[i*column + j];
            return r;
        }

        Matrix operator += (const Matrix& m)
        {
            if (row != m.row || column != m.column) {
                throw std::invalid_argument("Matrix addition: Matrix dimensions mismatch");
            }
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    (*this)(i, j) += m(i, j);
            return (*this);
        }

        Matrix operator - (const Matrix& m) const
        {
            if (row != m.row && column != m.column) {
                throw std::invalid_argument("Matrix subtraction: Matrix dimensions mismatch");
            }
            Matrix r(row, column);
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    r(i, j) = (*this)(i, j) - m(i, j);
            return r;
        }

        Matrix operator - (void) const
        {
            Matrix r(row, column);
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    r(i, j) = -(*this)(i, j);
            return r;
        }

        Matrix operator -= (const Matrix& m)
        {
            if (row != m.row || column != m.column) {
                throw std::invalid_argument("Matrix subtraction: Matrix dimensions mismatch");
            }
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    (*this)(i,j) -= m(i, j);
            return (*this);
        }

        Matrix operator * (const Matrix& m) const
        {
            if (column != m.row) {
                throw std::invalid_argument("Matrix multiplication: Matrix dimensions mismatch");
            }
            Matrix r(row, m.column);
            for (unsigned char i = 0; i < row; i++){
                int rRow = i*r.column; int thisRow = i*column;
                for (unsigned char k = 0; k < column; k++)
                {
                    float s = matrix[thisRow + k];
                    int mRow = k*m.column;
                    for (unsigned char j = 0; j < m.column; j++)
                        r.matrix[rRow + j] += (s * m.matrix[mRow + j]);
                }
            }
            return r;
        }
        //�������˷�
        Matrix<Scalar> operator * (const Vector3<Scalar>& vector) const
        {
            if (column != 3) {
                throw std::invalid_argument("Matrix multiplication: Matrix dimensions mismatch");
            }
            Matrix<Scalar> r(row, 1);
            r = (*this)*vector.toColVector();
            return r;
        }

        Matrix operator * (const float& n) const
        {
            Matrix r(row, column);
            for (unsigned char i = 0; i < row; i++){
                int rRow = i*r.column;int thisRow = i*column;
                for (unsigned char j = 0; j < column; j++)
                    r.matrix[rRow + j] = matrix[thisRow + j] * n;
            }
            return r;
        }

        Matrix operator * (const double& n) const
        {
            Matrix r(row, column);
            for (unsigned char i = 0; i < row; i++){
                int rRow = i*r.column;int thisRow = i*column;
                for (unsigned char j = 0; j < column; j++)
                    r.matrix[rRow + j] = matrix[thisRow + j] * n;
            }
            return r;
        }


        Matrix operator = (const Matrix& m)
        {
            // if (row == m.row && column == m.column)
            // {
            //     for (unsigned char i = 0; i < row; i++)
            //         for (unsigned char j = 0; j < column; j++)
            //             matrix[i * column + j] = m.matrix[i * column + j];
            //     return *this;
            // }
            // return *this;
            if (this != &m) {
                if (row != m.rows() || column != m.cols())//����ߴ粻һ�£������޸Ĵ�С
                {
                    delete[] matrix;
                    row = m.row;
                    column = m.column;
                    matrix = new Scalar[row * column];
                }
                std::copy(m.matrix, m.matrix + row * column, matrix);
            }
            return *this;
        }

        //������ȡ��������ķֿ����
        Matrix<Scalar>& operator = (const Block<Scalar>& m)
        {
            if (this == m.getParent())
                throw std::runtime_error("Block can not assign parent Matrix");
            if (row != m.rows() || column != m.cols())//����ߴ粻һ�£������޸Ĵ�С
            {
                delete[] matrix;
                row = m.rows();
                column = m.cols();
                matrix = new Scalar[row * column];
            }
            for (int i = 0; i < m.rows(); i++)
            {
                std::copy(&(m(i, 0)), &(m(i, 0)) + column, matrix+i*column);
            }
            return *this;
        }

        //������ȡvector
        Matrix<Scalar>& operator = (const Vector3<Scalar>& vector)
        {
            (*this) = vector.toColVector();
            return *this;
        }

        // float& operator () (unsigned char n)
        // {
        //     if (row == 1 || column == 1)
        //         return matrix[n];
        // }
        // float* operator [] (unsigned char n)
        // {
        //     if (n < row)
        //         return matrix + n * column;
        //     return matrix;
        // }
    public:
        void resize(int x, int y) //�ؽ�����
        {
            // delete[] mp;
            delete[] matrix;
            row = x;
            column = y;
            matrix = new Scalar[row * column];
            // mp = new Scalar * [row];
            // for (int i = 0; i < row; i++) {
            //     mp[i] = (matrix + i * column);
            // }
            setZero();
        }

        void setZero(void)//��������
        {
            std::fill(matrix, matrix + row * column, Scalar(0));
        }

        static Matrix Zero(int x, int y) //����ȫ0����
        {
            Matrix mat(x, y);
            mat.setZero();
            return mat;
        }

        Matrix setIdentity(void)//����λ��
        {
            for (unsigned char i = 0; i < row; i++)
            {
                for (unsigned char j = 0; j < column; j++)
                {
                    if (i == j)
                        matrix[i*column + j] = 1;
                    else
                        matrix[i*column + j] = 0;
                }
            }
            return *this;
        }

        static Matrix Identity(int rows, int cols)//������λ����
        {
            Matrix r(rows, cols);
            for (unsigned char i = 0; i < rows; i++) {
                for (unsigned char j = 0; j < cols; j++)
                {
                    if (i == j)
                        r.matrix[i*r.column + j] = 1;
                    else
                        r.matrix[i*r.column + j] = 0;
                }
            }
            return r;
        }


        void addIdentity(void)//�Լӵ�λ��
        {
            if (row == column)
            {
                for (unsigned char i = 0; i < row; i++)
                    matrix[i*column + i] += 1;
            }

        }

        void Transpose(void)//����ת��
        {
            if (row == column)
            {
                Scalar t;
                for (unsigned char i = 0; i < row - 1; i++)
                    for (unsigned char j = i + 1; j < column; j++)
                    {
                        t = matrix[i*column + j];
                        matrix[i*column + j] = matrix[j*column + i];
                        matrix[j*column + i] = t;
                    }
            }
            else
            {
                Matrix t(column, row);
                for (unsigned char i = 0; i < row; i++)
                    for (unsigned char j = 0; j < column; j++)
                        t.matrix[j*t.column + i] = matrix[i*column + j];
                (*this) = t;
            }

        }

        Matrix transpose(void)//����ת�þ���
        {
            Matrix r(column, row);
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    r.matrix[j*r.column + i] = matrix[i*column + j];
            return r;

        }

        void Inverse()//�������
        {
            if (row == column)
            {
                Matrix t(row, column);
                Matrix Cofactor(row - 1, row - 1);

                for (unsigned char i = 0; i < row; i++)
                    for (unsigned char j = 0; j < column; j++)
                    {
                        Cofactor = adjoint(*this, i, j);
                        t.matrix[j*t.column + i] = (1 - 2 * ((i + j) % 2)) * determinant(Cofactor);
                    }
                t = t * (1 / determinant(*this));

                for (unsigned char i = 0; i < row; i++)
                    for (unsigned char j = 0; j < column; j++)
                        matrix[i*column + j] = t.matrix[i*t.column + j];


            }
        }

        Matrix inverse()//���������
        {
            if (row == column)
            {
                float a;
                Matrix r(row, column);
                Matrix Cofactor(row - 1, row - 1);

                for (unsigned char i = 0; i < row; i++)
                    for (unsigned char j = 0; j < column; j++)
                    {
                        Cofactor = adjoint(*this, i, j);
                        r.matrix[j*r.column + i] = (1 - 2 * ((i + j) % 2)) * determinant(Cofactor);
                    }
                a = determinant(*this);
                if (a)
                    r = r * (1 / a);

                return r;
            }
            Matrix r(row, column);
            r.setZero();
            return r;
        }


        void set_Block_Matrix(Matrix& m, unsigned char r, unsigned char c)//���ֿ����m�����Ӧλ�ã�m(0,0)��Ӧr��c
        {
            if (m.row + r <= row || m.column + c <= column)
            {
                for (unsigned char i = 0; i < m.row; i++)
                    for (unsigned char j = 0; j < m.column; j++)
                        matrix[(r + i) * column + j + c] = m.matrix[i * m.column + j];
            }
        }

        // ���߽���ĳ������汾
        //P,QΪ�ֿ��С
        template<int P, int Q>
        Block<Scalar> block(int startRow, int startCol) {
            static_assert(P > 0 && Q > 0, "Block size must be positive");
            return block(startRow, startCol, P, Q);
        }

        Block<Scalar> block(int startRow, int startCol, int blockRows, int blockCols)//ȥ����ʼ��λ��r��c��h��w����Ϊ�ֿ����
        {
            //��֤���������ȷ��
            if (startRow < 0 || startCol < 0 ||
                blockRows <= 0 || blockCols <= 0 ||
                startRow + blockRows > row ||
                startCol + blockCols > column)
            {
                throw std::out_of_range(
                    "Matrix block indices out of range: [" +
                    std::to_string(startRow) + "," + std::to_string(startCol) + "]+" +
                    std::to_string(blockRows) + "x" + std::to_string(blockCols) +
                    " in " + std::to_string(row) + "x" + std::to_string(column) + " matrix"
                );
            }
            //Matrix m(blockRows, blockCols);
            //for (unsigned char i = 0; i < blockRows; i++)
            //    for (unsigned char j = 0; j < blockCols; j++)
            //        m.matrix[i * blockCols + j] = matrix[(startRow + i) * column + j + startCol];
            Block<Scalar> m(*this, startRow, startCol, blockRows, blockCols);
            return m;
        }


        float trace()//�����ļ�
        {
            if (row != column)
                throw std::invalid_argument("Matrix trace: Matrix row not equal column");

            float t = 0;
            for (unsigned char i = 0; i < row; i++)
                t += (*this)(i, i);
            return t;
        }



    public:
        static Matrix adjoint(Matrix& m, unsigned char r, unsigned char c)//��������
        {
            if (m.row == m.column && r < m.row && c < m.column)
            {
                Matrix cofactor(m.row - 1, m.row - 1);

                unsigned char x, y;
                for (unsigned char i = 0; i < m.row - 1; i++)
                    for (unsigned char j = 0; j < m.row - 1; j++)
                    {
                        if (i >= r && r != m.row)
                            x = i + 1;
                        else
                            x = i;

                        if (j >= c && c != m.row)
                            y = j + 1;
                        else
                            y = j;

                        cofactor(i, j) = m(x, y);
                    }
                return cofactor;
            }
            Matrix la(m.row - 1, m.row - 1);
            la.setZero();
            return la;
        }



        static float determinant(Matrix& m)//������ʽ
        {

            if (m.row == m.column)
            {
                float result = 0, cofactor = 0;
                short sign = 1;
                Matrix Cofactor(m.row - 1, m.row - 1);
                if (m.row == 1)
                    result = m.matrix[0];
                else if (m.row == 2)
                {
                    result = m.matrix[0] * m.matrix[3] - m.matrix[1] * m.matrix[2];
                }
                else
                {
                    for (unsigned char i = 0; i < m.row; i++)
                    {
                        cofactor = m.matrix[i * m.column];
                        if (cofactor != 0)
                        {
                            Cofactor = adjoint(m, i, 0);
                            result += sign * determinant(Cofactor) * cofactor;
                        }
                        sign *= -1;
                    }
                }
                return result;
            }
            float result = 0;
            return result;
        }



    };

    // ���� Eigen ���ı���
    using MatrixXf = Matrix<float>;
    using MatrixXd = Matrix<double>;
}
#endif /* IMU_MATRIX_H_ */
