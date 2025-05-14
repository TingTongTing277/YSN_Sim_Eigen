/*
 * Vector.h
 *
 *  Created on: 2023年4月9日
 *      Author: fuzic
 */
#pragma once
#ifndef IMU_VECTOR_H_
#define IMU_VECTOR_H_

#include "Matrix.h"

namespace Sim_Eigen {
    template <typename Scalar>

    class Vector3
    {
        
   private:
       Scalar v[3];
    public:
        Vector3() {
            // 通过 this-> 访问基类成员（避免模板依赖问题）
            this->v[0] = 0;
            this->v[1] = 0;
            this->v[2] = 0;
        }
        Vector3(Scalar x, Scalar y, Scalar z) {
            // 通过 this-> 访问基类成员（避免模板依赖问题）
            this->v[0] = x;
            this->v[1] = y;
            this->v[2] = z;
        }
        Vector3(const Vector3& other) {
            this->v[0] = other(0);
            this->v[1] = other(1);
            this->v[2] = other(2);
        }

        ~Vector3() {};

        // Eigen兼容接口
        static constexpr int SizeAtCompileTime = 3;
        using ScalarType = Scalar;

        template<typename, typename> friend class CommaInitializer;

        template<typename Derived, typename T2>
        class CommaInitializer {
            Derived& vector_;
            int index_;
        public:
            CommaInitializer(Derived& vec, T2 first_val) : vector_(vec), index_(0) {
                vector_(index_++) = static_cast<typename Derived::ScalarType>(first_val);
            }

            CommaInitializer& operator,(T2 value) {
                if (index_ >= 3) {
                    throw std::out_of_range("Too many elements for Vector3 initialization");
                }
                vector_(index_++) = static_cast<typename Derived::ScalarType>(value);
                return *this;
            }

            ~CommaInitializer() noexcept(false) {
                if (index_ != 3) {
                    throw std::invalid_argument(
                        "Vector3 requires exactly 3 elements, got " +
                        std::to_string(index_));
                }
            }
        };

    public:

        // 元素访问
        inline Scalar& operator()(int index) {
            if (index < 0 || index >= 3)
                throw std::out_of_range("Vector3 index out of range");
            return v[index];
        }

        inline const Scalar& operator()(int index) const {
            if (index < 0 || index >= 3)
                throw std::out_of_range("Vector3 index out of range");
            return v[index];
        }

        // 元素访问
        inline Scalar& operator[](int index) {
            if (index < 0 || index >= 3)
                throw std::out_of_range("Vector3 index out of range");
            return v[index];
        }

        inline const Scalar& operator[](int index) const {
            if (index < 0 || index >= 3)
                throw std::out_of_range("Vector3 index out of range");
            return v[index];
        }

        Scalar& x() { return v[0]; }
        Scalar& y() { return v[1]; }
        Scalar& z() { return v[2]; }
        const Scalar& x() const { return v[0]; }
        const Scalar& y() const { return v[1]; }
        const Scalar& z() const { return v[2]; }

        Vector3 operator()(Scalar x, Scalar y, Scalar z)
        {
            v[0] = x;
            v[1] = y;
            v[2] = z;
            return*this;
        }

        Vector3 operator = (const Matrix<Scalar>& m)//向量赋值
        {
            if(m.rows() != 3 || m.cols() != 1)
                throw std::runtime_error("Matrix size is not equal vector size");
            this->v[0] = m(0,0);
            this->v[1] = m(1,0);
            this->v[2] = m(2,0);
            return *this;
        }
        
        //显式定义拷贝赋值运算符
        Vector3& operator = (const Vector3& vector)//向量赋值
        {
            (*this)(0) = vector(0); (*this)(1) = vector(1); (*this)(2) = vector(2);// 调用基类深拷贝
            return *this;
        }
        // 显式定义移动赋值运算符（可选）
        Vector3& operator = (Vector3&& vector) noexcept {
            (*this)(0) = vector(0); (*this)(1) = vector(1); (*this)(2) = vector(2);
            return *this;
        }
        
        Vector3 operator + (const Vector3& vector) const//向量相加
        {
            Vector3 r;
            r.v[0] = vector.v[0] + v[0];
            r.v[1] = vector.v[1] + v[1];
            r.v[2] = vector.v[2] + v[2];
            return r;
        }

        Vector3 operator +=(const Vector3& vector)//向量自加
        {
            this->v[0] += vector.v[0];
            this->v[1] += vector.v[1];
            this->v[2] += vector.v[2];
            return *this;
        }


        Vector3 operator - (const Vector3& vector) const//向量相减
        {
            Vector3 r;
            r.v[0] = this->v[0] - vector.v[0];
            r.v[1] = this->v[1] - vector.v[1];
            r.v[2] = this->v[2] - vector.v[2];
            return r;
        }

        Vector3 operator - () const//向量相反
        {
            Vector3 r;
            r.v[0] = -(this->v[0]);
            r.v[1] = -(this->v[1]);
            r.v[2] = -(this->v[2]);
            return r;

        }
        
        Vector3 operator -= (const Vector3& vector)//向量自减
        {
            this->v[0] -= vector.v[0];
            this->v[1] -= vector.v[1];
            this->v[2] -= vector.v[2];
            return *this;
        }

        /*
        Scalar  operator * (const Vector3& vector) const//点乘
        {
            Scalar r = v[0] * vector(0) + v[1] * vector(1) + v[2] * vector(2);
            return r;
        }*/
        Vector3 operator * (const Scalar& n) const //向量数乘
        {
            Vector3 r((*this)(0) * n, (*this)(1) * n, (*this)(2) * n);
            return r;
        }

        bool operator == (const Vector3& other) const
        {
            if((*this)(0)==other(0)&&(*this)(1)==other(1)&&(*this)(2)==other(2))
                return true;
            else
                return false;
        }
        
    public:
        Scalar norm(void)//取模
        {
            Scalar r = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return r;
        }

        Scalar squaredNorm() const {
            return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        }

        void normalize(void)//标准化
        {
            Scalar t = this->norm();
            if (t != 0)
            {
                (*this)(0) = (*this)(0) / t;
                (*this)(1) = (*this)(1) / t;
                (*this)(2) = (*this)(2) / t;
            }
            else
                throw std::runtime_error("Cannot normalize zero vector");
        }

        Vector3 normalized(void)//返回标准化（不改变原值）
        {
            Vector3 r;
            Scalar t = this->norm();
            if (t != 0)
            {
                r(0) = (*this)(0) / t;
                r(1) = (*this)(1) / t;
                r(2) = (*this)(2) / t;
            }
            else
                throw std::runtime_error("Cannot normalize zero vector");

            return r;
        }

        Scalar dot(const Vector3& vector)//点乘
        {
            Scalar r = (*this)(0) * vector(0) + (*this)(1) * vector(1) + (*this)(2) * vector(2);
            return r;
        }

        Vector3 cross(const Vector3& vector)//叉乘
        {
            Vector3 r;
            r(0) = (*this)(1) * vector(2) - (*this)(2) * vector(1);
            r(1) = (*this)(2) * vector(0) - (*this)(0) * vector(2);
            r(2) = (*this)(0) * vector(1) - (*this)(1) * vector(0);
            return r;
        }

        Vector3 get_cwiseInverse()//各个位置上的倒数
        {
            Vector3 r;
            r(0) = 1 / (*this)(0);
            r(1) = 1 / (*this)(1);
            r(2) = 1 / (*this)(2);

            return r;

        }

        Vector3 cwiseProduct(const Vector3& vector)//各个位置上的相乘
        {
            Vector3 r;
            r(0) = (*this)(0) * vector(0);
            r(1) = (*this)(1) * vector(1);
            r(2) = (*this)(2) * vector(2);
            return r;
        }

        Vector3 cwiseQuotient(const Vector3& other) const {//各个位置上的相除
            return Vector3((*this)(0) / other(0), (*this)(1) / other(1), (*this)(2) / other(2));
        }

        Matrix<Scalar> get_Asdiagonal(void)//向量变对角阵
        {
            Matrix<Scalar> r(3, 3);
            r(0, 0) = (*this)(0);
            r(1, 1) = (*this)(1);
            r(2, 2) = (*this)(2);
            return r;
        }

        Matrix<Scalar> get_Antisymmetric_Matrix(void)//求反对称矩阵
        {
            Matrix<Scalar> r(3, 3);
            r(0, 0) = 0, r(0, 1) = -(*this)(2), r(0, 2) = (*this)(1);
            r(1, 0) = (*this)(2), r(1, 1) = 0, r(1, 2) = -(*this)(0);
            r(2, 0) = -(*this)(1), r(2, 1) = (*this)(0), r(2, 2) = 0;
            return r;
        }



        void setZero(void)//清0
        {
            (*this)(0) = 0;
            (*this)(1) = 0;
            (*this)(2) = 0;
        }
        void setOnes(void)//变1（全元素）
        {
            (*this)(0) = 1;
            (*this)(1) = 1;
            (*this)(2) = 1;
        }


        Matrix<Scalar> toRowVector(void) const//转矩阵行向量
        {
            Matrix<Scalar> r(1, 3);
            r(0, 0) = (*this)(0), r(0, 1) = (*this)(1), r(0, 2) = (*this)(2);
            return r;
        }

        Matrix<Scalar> toColVector(void) const//转矩阵列向量
        {
            Matrix<Scalar> r(3, 1);
            r(0, 0) = (*this)(0), r(1, 0) = (*this)(1), r(2, 0) = (*this)(2);
            return r;

        }

        // 流式初始化（需要Matrix.h已实现逗号初始化）
        CommaInitializer<Vector3<Scalar>, Scalar> operator<<(Scalar value) {
            return CommaInitializer<Vector3<Scalar>, Scalar>(*this, value);
        }

        // 添加友元声明
        friend std::ostream& operator<<(std::ostream& os, const Vector3<Scalar>& vector) {
            // 保存原有格式状态
            std::ios_base::fmtflags origFlags = os.flags();
            std::streamsize origPrecision = os.precision();
            std::streamsize origWidth = os.width();

            // 设置输出格式（示例：固定小数点，4位精度）
            os << std::fixed << std::setprecision(4);

            for (int i = 0; i < 3; ++i) {
                os << vector(i);
                os << std::endl;  // 行末换行
            }

            // 恢复原有格式
            os.flags(origFlags);
            os.precision(origPrecision);
            os.width(origWidth);
            return os;
        }

        // 特殊操作
        static Vector3 Zero() { return Vector3(0, 0, 0); }
        static Vector3 Ones() { return Vector3(1, 1, 1); }
        static Vector3 UnitX() { return Vector3(1, 0, 0); }
        static Vector3 UnitY() { return Vector3(0, 1, 0); }
        static Vector3 UnitZ() { return Vector3(0, 0, 1); }

    };


    typedef Vector3<float> Vector3f;
    typedef Vector3<double> Vector3d;


}
#endif /* IMU_VECTOR_H_ */
