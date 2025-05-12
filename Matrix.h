/*
 * Matrix.h
 *
 *  Created on: 2023年4月9日
 *      Author: fuzic
 */
#pragma once
#ifndef MATH_MATRIX_H_
#define MATH_MATRIX_H_

//#include "Math_Header.h"
#include <stdexcept>  // 添加异常处理
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

        // 逗号初始化器实现
        template<typename, typename> friend class CommaInitializer;

        template<typename Derived, typename Scalar2>
        class CommaInitializer {
            Derived& matrix_;
            int index_;
            int expectedSize_;
        public:
            CommaInitializer(Derived& mat, Scalar2 first_val) : matrix_(mat), index_(0) {
                expectedSize_ = mat.rows() * mat.cols();
                *this, first_val;  // 处理第一个值
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
        Matrix(int x, int y) :row(x), column(y)//矩阵定义
        {
            if (x > 0 && y > 0) {
                matrix = new Scalar[row * column];
                setZero();
            }
            else {
                matrix = nullptr;
            }
        }
        Matrix(Matrix& other) :row(other.row), column(other.column)//拷贝构造函数
        {
            matrix = new Scalar[row * column];
            std::copy(other.matrix, other.matrix + row * column, matrix);
        }

        // 移动构造函数
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
        // Eigen 风格接口
        int rows() const { return row; }
        int cols() const { return column; }
        Scalar* data() { return matrix; }
        const Scalar* data() const { return matrix; }

        // 添加友元声明
        friend std::ostream& operator<<(std::ostream& os, const Matrix<Scalar>& mat) {
            // 保存原有格式状态
            std::ios_base::fmtflags origFlags = os.flags();
            std::streamsize origPrecision = os.precision();
            std::streamsize origWidth = os.width();

            // 设置输出格式（示例：固定小数点，4位精度）
            os << std::fixed << std::setprecision(4);

            for (int i = 0; i < mat.rows(); ++i) {
                for (int j = 0; j < mat.cols(); ++j) {
                    os << mat(i, j);
                    if (j != mat.cols() - 1) os << "\t";  // 元素间用空格分隔
                }
                if (i != mat.rows() - 1) os << std::endl;  // 行末换行
            }
            os << std::endl; // 行末换行

            // 恢复原有格式
            os.flags(origFlags);
            os.precision(origPrecision);
            os.width(origWidth);
            return os;
        }

        // 运算符重载
        // 流式初始化运算符
        CommaInitializer<Matrix<Scalar>, Scalar> operator<<(Scalar value) {
            if (!matrix) throw std::runtime_error("Matrix not initialized");
            return CommaInitializer<Matrix<Scalar>, Scalar>(*this, value);
        }

        //括号重载
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
        //对向量乘法
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
                if (row != m.rows() || column != m.cols())//矩阵尺寸不一致，重新修改大小
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

        //用于提取其他矩阵的分块矩阵
        Matrix<Scalar>& operator = (const Block<Scalar>& m)
        {
            if (this == m.getParent())
                throw std::runtime_error("Block can not assign parent Matrix");
            if (row != m.rows() || column != m.cols())//矩阵尺寸不一致，重新修改大小
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

        //用于提取vector
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
        void resize(int x, int y) //重建矩阵
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

        void setZero(void)//矩阵清零
        {
            std::fill(matrix, matrix + row * column, Scalar(0));
        }

        static Matrix Zero(int x, int y) //建立全0矩阵
        {
            Matrix mat(x, y);
            mat.setZero();
            return mat;
        }

        Matrix setIdentity(void)//矩阵单位化
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

        static Matrix Identity(int rows, int cols)//建立单位矩阵
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


        void addIdentity(void)//自加单位阵
        {
            if (row == column)
            {
                for (unsigned char i = 0; i < row; i++)
                    matrix[i*column + i] += 1;
            }

        }

        void Transpose(void)//矩阵转置
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

        Matrix transpose(void)//返回转置矩阵
        {
            Matrix r(column, row);
            for (unsigned char i = 0; i < row; i++)
                for (unsigned char j = 0; j < column; j++)
                    r.matrix[j*r.column + i] = matrix[i*column + j];
            return r;

        }

        void Inverse()//求逆矩阵
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

        Matrix inverse()//返回逆矩阵
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


        void set_Block_Matrix(Matrix& m, unsigned char r, unsigned char c)//将分块矩阵m放入对应位置，m(0,0)对应r，c
        {
            if (m.row + r <= row || m.column + c <= column)
            {
                for (unsigned char i = 0; i < m.row; i++)
                    for (unsigned char j = 0; j < m.column; j++)
                        matrix[(r + i) * column + j + c] = m.matrix[i * m.column + j];
            }
        }

        // 带边界检查的常量化版本
        //P,Q为分块大小
        template<int P, int Q>
        Block<Scalar> block(int startRow, int startCol) {
            static_assert(P > 0 && Q > 0, "Block size must be positive");
            return block(startRow, startCol, P, Q);
        }

        Block<Scalar> block(int startRow, int startCol, int blockRows, int blockCols)//去除起始点位于r、c的h行w列作为分块矩阵
        {
            //验证输入参数正确性
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


        float trace()//求矩阵的迹
        {
            if (row != column)
                throw std::invalid_argument("Matrix trace: Matrix row not equal column");

            float t = 0;
            for (unsigned char i = 0; i < row; i++)
                t += (*this)(i, i);
            return t;
        }



    public:
        static Matrix adjoint(Matrix& m, unsigned char r, unsigned char c)//求伴随矩阵
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



        static float determinant(Matrix& m)//求行列式
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

    // 定义 Eigen 风格的别名
    using MatrixXf = Matrix<float>;
    using MatrixXd = Matrix<double>;
}
#endif /* IMU_MATRIX_H_ */
