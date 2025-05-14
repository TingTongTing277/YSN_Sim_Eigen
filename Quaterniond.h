/*
 * Quaterniond.h
 *
 *  Created on: 2023年4月9日
 *      Author: fuzic
 */
#pragma once
#ifndef IMU_QUATERNIOND_H_
#define IMU_QUATERNIOND_H_


#include "Vector.h"

namespace Sim_Eigen {
    template<typename Scalar>
    class Quaternion
    {

    private:
        Scalar w_;
        Scalar x_;
        Scalar y_;
        Scalar z_;

    public:
        Quaternion() :w_(0), x_(0), y_(0), z_(0) {}
        Quaternion(Scalar wi, Scalar xi, Scalar yi, Scalar zi) :w_(wi), x_(xi), y_(yi), z_(zi) {}
        Quaternion(const Quaternion& q)
        {
            w_ = q.w();
            x_ = q.x();
            y_ = q.y();
            z_ = q.z();
        }


    public:
        Scalar& x() { return this->x_; }
        Scalar& y() { return this->y_; }
        Scalar& z() { return this->z_; }
        Scalar& w() { return this->w_; }
        const Scalar& x() const { return this->x_; }
        const Scalar& y() const { return this->y_; }
        const Scalar& z() const { return this->z_; }
        const Scalar& w() const { return this->w_; }

        Quaternion operator =(const Quaternion& q)
        {
            w_ = q.w();
            x_ = q.x();
            y_ = q.y();
            z_ = q.z();
            return *this;
        }

        Quaternion operator + (const Quaternion& q) const
        {
            Quaternion r;
            r.w() = w_ + q.w();
            r.x() = x_ + q.x();
            r.y() = y_ + q.y();
            r.z() = z_ + q.z();
            return r;
        }

        Quaternion operator - (const Quaternion& q) const
        {
            Quaternion r;
            r.w() = w_ - q.w();
            r.x() = x_ - q.x();
            r.y() = y_ - q.y();
            r.z() = z_ - q.z();
            return r;
        }
        Quaternion operator - () const
        {
            Quaternion r;
            r.w() = -w_;
            r.x() = -x_;
            r.y() = -y_;
            r.z() = -z_;
            return r;
        }




        Quaternion operator * (const Quaternion& q) const
        {
            Quaternion r;

            r.w() = w_ * q.w() - x_ * q.x() - y_ * q.y() - z_ * q.z();
            r.x() = w_ * q.x() + x_ * q.w() + y_ * q.z() - z_ * q.y();
            r.y() = w_ * q.y() + y_ * q.w() + z_ * q.x() - x_ * q.z();
            r.z() = w_ * q.z() + z_ * q.w() + x_ * q.y() - y_ * q.x();

            return r;
        }


        Quaternion operator * (const float n) const
        {
            Quaternion r;
            r.w() = n * w_;
            r.x() = n * x_;
            r.y() = n * y_;
            r.z() = n * z_;
            return r;
        }



    public:
        float norm(void) const
        {
            float r;
            r = sqrt(w_ * w_ + x_ * x_ + y_ * y_ + z_ * z_);
            return r;
        }

        void normalize(void)
        {
            float t = this->norm();
            if (t != 0)
            {
                w_ = w_ / t;
                x_ = x_ / t;
                y_ = y_ / t;
                z_ = z_ / t;
            }
            else
                throw std::runtime_error("Cannot normalize zero Quaternion");
        }

        Quaternion normalized(void)
        {
            float t = this->norm();
            if (t != 0)
            {
                Quaternion q;
                q.w() = w_ / t;
                q.x() = x_ / t;
                q.y() = y_ / t;
                q.z() = z_ / t;
                return q;
            }
            else
                throw std::runtime_error("Cannot normalize zero Quaternion");
        }

        Scalar dot(Quaternion other)
        {
            Scalar result;
            result = w_*other.w() + x_*other.x() + y_*other.y() + z_*other.z();
            return result;
        }

        Quaternion conjugate(void)
        {
            x_ = -x_;
            y_ = -y_;
            z_ = -z_;
            return (*this);
        }

        Quaternion get_conjugate(void)
        {
            Quaternion q;
            q.w() = w_;
            q.x() = -x_;
            q.y() = -y_;
            q.z() = -z_;
            return q;
        }
        Quaternion inverse(void)
        {
            float t = this->norm();
            if (t)
                *this = this->get_conjugate() * (1 / (t * t));
            else
                this->w_ = 1;
            return (*this);
        }

        Quaternion get_inverse(void)
        {
            float t = this->norm();
            Quaternion q;
            if (t)
                q = this->get_conjugate() * (1 / (t * t));
            else
                q.w_ = 1;
            return q;

        }

        Vector3<Scalar> vec(void)
        {
            Vector3<Scalar> r(x_, y_, z_);
            return r;
        }


        float get_cos_half(void)
        {
            return  w_ / (this->norm());
        }

        float get_sin_half(void)
        {
            return sqrt(x_ * x_ + y_ * y_ + z_ * z_) / (this->norm());
        }




    };

    // 定义 Eigen 风格的别名
    using Quaternionf = Quaternion<float>;
    using Quaterniond = Quaternion<double>;

}
#endif /* IMU_QUATERNIOND_H_ */
