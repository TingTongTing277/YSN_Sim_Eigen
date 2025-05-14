/*
 * Vector.h
 *
 *  Created on: 2023��4��9��
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
            // ͨ�� this-> ���ʻ����Ա������ģ���������⣩
            this->v[0] = 0;
            this->v[1] = 0;
            this->v[2] = 0;
        }
        Vector3(Scalar x, Scalar y, Scalar z) {
            // ͨ�� this-> ���ʻ����Ա������ģ���������⣩
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

        // Eigen���ݽӿ�
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

        // Ԫ�ط���
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

        // Ԫ�ط���
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

        Vector3 operator = (const Matrix<Scalar>& m)//������ֵ
        {
            if(m.rows() != 3 || m.cols() != 1)
                throw std::runtime_error("Matrix size is not equal vector size");
            this->v[0] = m(0,0);
            this->v[1] = m(1,0);
            this->v[2] = m(2,0);
            return *this;
        }
        
        //��ʽ���忽����ֵ�����
        Vector3& operator = (const Vector3& vector)//������ֵ
        {
            (*this)(0) = vector(0); (*this)(1) = vector(1); (*this)(2) = vector(2);// ���û������
            return *this;
        }
        // ��ʽ�����ƶ���ֵ���������ѡ��
        Vector3& operator = (Vector3&& vector) noexcept {
            (*this)(0) = vector(0); (*this)(1) = vector(1); (*this)(2) = vector(2);
            return *this;
        }
        
        Vector3 operator + (const Vector3& vector) const//�������
        {
            Vector3 r;
            r.v[0] = vector.v[0] + v[0];
            r.v[1] = vector.v[1] + v[1];
            r.v[2] = vector.v[2] + v[2];
            return r;
        }

        Vector3 operator +=(const Vector3& vector)//�����Լ�
        {
            this->v[0] += vector.v[0];
            this->v[1] += vector.v[1];
            this->v[2] += vector.v[2];
            return *this;
        }


        Vector3 operator - (const Vector3& vector) const//�������
        {
            Vector3 r;
            r.v[0] = this->v[0] - vector.v[0];
            r.v[1] = this->v[1] - vector.v[1];
            r.v[2] = this->v[2] - vector.v[2];
            return r;
        }

        Vector3 operator - () const//�����෴
        {
            Vector3 r;
            r.v[0] = -(this->v[0]);
            r.v[1] = -(this->v[1]);
            r.v[2] = -(this->v[2]);
            return r;

        }
        
        Vector3 operator -= (const Vector3& vector)//�����Լ�
        {
            this->v[0] -= vector.v[0];
            this->v[1] -= vector.v[1];
            this->v[2] -= vector.v[2];
            return *this;
        }

        /*
        Scalar  operator * (const Vector3& vector) const//���
        {
            Scalar r = v[0] * vector(0) + v[1] * vector(1) + v[2] * vector(2);
            return r;
        }*/
        Vector3 operator * (const Scalar& n) const //��������
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
        Scalar norm(void)//ȡģ
        {
            Scalar r = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return r;
        }

        Scalar squaredNorm() const {
            return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        }

        void normalize(void)//��׼��
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

        Vector3 normalized(void)//���ر�׼�������ı�ԭֵ��
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

        Scalar dot(const Vector3& vector)//���
        {
            Scalar r = (*this)(0) * vector(0) + (*this)(1) * vector(1) + (*this)(2) * vector(2);
            return r;
        }

        Vector3 cross(const Vector3& vector)//���
        {
            Vector3 r;
            r(0) = (*this)(1) * vector(2) - (*this)(2) * vector(1);
            r(1) = (*this)(2) * vector(0) - (*this)(0) * vector(2);
            r(2) = (*this)(0) * vector(1) - (*this)(1) * vector(0);
            return r;
        }

        Vector3 get_cwiseInverse()//����λ���ϵĵ���
        {
            Vector3 r;
            r(0) = 1 / (*this)(0);
            r(1) = 1 / (*this)(1);
            r(2) = 1 / (*this)(2);

            return r;

        }

        Vector3 cwiseProduct(const Vector3& vector)//����λ���ϵ����
        {
            Vector3 r;
            r(0) = (*this)(0) * vector(0);
            r(1) = (*this)(1) * vector(1);
            r(2) = (*this)(2) * vector(2);
            return r;
        }

        Vector3 cwiseQuotient(const Vector3& other) const {//����λ���ϵ����
            return Vector3((*this)(0) / other(0), (*this)(1) / other(1), (*this)(2) / other(2));
        }

        Matrix<Scalar> get_Asdiagonal(void)//������Խ���
        {
            Matrix<Scalar> r(3, 3);
            r(0, 0) = (*this)(0);
            r(1, 1) = (*this)(1);
            r(2, 2) = (*this)(2);
            return r;
        }

        Matrix<Scalar> get_Antisymmetric_Matrix(void)//�󷴶Գƾ���
        {
            Matrix<Scalar> r(3, 3);
            r(0, 0) = 0, r(0, 1) = -(*this)(2), r(0, 2) = (*this)(1);
            r(1, 0) = (*this)(2), r(1, 1) = 0, r(1, 2) = -(*this)(0);
            r(2, 0) = -(*this)(1), r(2, 1) = (*this)(0), r(2, 2) = 0;
            return r;
        }



        void setZero(void)//��0
        {
            (*this)(0) = 0;
            (*this)(1) = 0;
            (*this)(2) = 0;
        }
        void setOnes(void)//��1��ȫԪ�أ�
        {
            (*this)(0) = 1;
            (*this)(1) = 1;
            (*this)(2) = 1;
        }


        Matrix<Scalar> toRowVector(void) const//ת����������
        {
            Matrix<Scalar> r(1, 3);
            r(0, 0) = (*this)(0), r(0, 1) = (*this)(1), r(0, 2) = (*this)(2);
            return r;
        }

        Matrix<Scalar> toColVector(void) const//ת����������
        {
            Matrix<Scalar> r(3, 1);
            r(0, 0) = (*this)(0), r(1, 0) = (*this)(1), r(2, 0) = (*this)(2);
            return r;

        }

        // ��ʽ��ʼ������ҪMatrix.h��ʵ�ֶ��ų�ʼ����
        CommaInitializer<Vector3<Scalar>, Scalar> operator<<(Scalar value) {
            return CommaInitializer<Vector3<Scalar>, Scalar>(*this, value);
        }

        // �����Ԫ����
        friend std::ostream& operator<<(std::ostream& os, const Vector3<Scalar>& vector) {
            // ����ԭ�и�ʽ״̬
            std::ios_base::fmtflags origFlags = os.flags();
            std::streamsize origPrecision = os.precision();
            std::streamsize origWidth = os.width();

            // ���������ʽ��ʾ�����̶�С���㣬4λ���ȣ�
            os << std::fixed << std::setprecision(4);

            for (int i = 0; i < 3; ++i) {
                os << vector(i);
                os << std::endl;  // ��ĩ����
            }

            // �ָ�ԭ�и�ʽ
            os.flags(origFlags);
            os.precision(origPrecision);
            os.width(origWidth);
            return os;
        }

        // �������
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
