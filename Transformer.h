/*
 * Transformer.h
 *
 *  Created on: 2023年4月9日
 *      Author: fuzic
 */
#pragma once
#ifndef MATH_TRANSFORMER_H_
#define MATH_TRANSFORMER_H_

#ifdef Eigen
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"
using namespace Eigen;
#else
#include "Quaterniond.h"
using namespace Sim_Eigen;
#endif // Eigen



class Rotation
{

#define M_PI 3.14159265358979323846

public:
    static Vector3f get_cwiseProduct(const Vector3f& v,const Vector3f& vector)//向量各个位置上的相乘
    {
        Vector3f r;
        r(0)=v(0)*vector(0);
        r(1)=v(1)*vector(1);
        r(2)=v(2)*vector(2);
        return r;
    }
    static MatrixXf get_Asdiagonal(const Vector3f& v)//向量变对角阵
    {
        MatrixXf r(3, 3);
        r.setZero();
        r(0,0)=v(0);
        r(1,1)=v(1);
        r(2,2)=v(2);
        return r;
    }
    static Vector3f matrix2euler(MatrixXf& dcm)//方向余弦阵转欧拉角
    {
        Vector3f euler;

        euler(1) = atan(-dcm(2,0) / sqrt(dcm(2,1) * dcm(2,1) + dcm(2,2) * dcm(2,2)));

        if (dcm(2,0) <= -0.999) {
            euler(0) = atan2(dcm(2,1), dcm(2,2));
            euler(2) = atan2((dcm(1,2) - dcm(0,1)), (dcm(0,2) + dcm(1,1)));
        }
        else if (dcm(2,0) >= 0.999) {
            euler(0) = atan2(dcm(2,1), dcm(2,2));
            euler(2) = (3.14159265358979323846 + atan2((dcm(1,2) + dcm(0,1)), (dcm(0,2) - dcm(1,1))));
        }
        else {
            euler(0) = atan2(dcm(2,1), dcm(2,2));
            euler(2) = atan2(dcm(1,0), dcm(0,0));
        }

        // heading 0~2PI
        if (euler(2) < 0) {
            euler(2) = 2 * 3.14159265358979323846 + euler(2);
        }

        return euler;
    }

    static MatrixXf euler2matrix(Vector3f& v)//欧拉角转方向余弦阵
    {
        MatrixXf r(3, 3);
        r(0,0) = cos(v(2)) * cos(v(1)), r(0,1) = -sin(v(2)) * cos(v(0)) + cos(v(2)) * sin(v(1)) * sin(v(0)), r(0,2) = sin(v(2)) * sin(v(0)) + cos(v(2)) * sin(v(1)) * cos(v(0));
        r(1,0) = sin(v(2)) * cos(v(1)), r(1,1) = cos(v(2)) * cos(v(0)) + sin(v(1)) * sin(v(0))*sin(v(2)), r(1,2) = -cos(v(2)) * sin(v(0)) + sin(v(1)) * sin(v(2)) * cos(v(0));
        r(2,0) = -sin(v(1)), r(2,1) = cos(v(1)) * sin(v(0)), r(2,2) = cos(v(1)) * cos(v(0));
        return r;
    }


    static Quaternionf  euler2quaternion(Vector3f& v)//欧拉角转四元数
    {
        Quaternionf r;

        r.w() = cos(v(2) / 2) * cos(v(1) / 2) * cos(v(0) / 2) + sin(v(2) / 2) * sin(v(1) / 2) * sin(v(0) / 2);
        r.z() = sin(v(2) / 2) * cos(v(1) / 2) * cos(v(0) / 2) - cos(v(2) / 2) * sin(v(1) / 2) * sin(v(0) / 2);
        r.y() = cos(v(2) / 2) * sin(v(1) / 2) * cos(v(0) / 2) + sin(v(2) / 2) * cos(v(1) / 2) * sin(v(0) / 2);
        r.x() = cos(v(2) / 2) * cos(v(1) / 2) * sin(v(0) / 2) - sin(v(2) / 2) * sin(v(1) / 2) * cos(v(0) / 2);

        return r;
    }


    static MatrixXf quaternion2matrix(Quaternionf& q)//欧拉角转方向余弦阵
    {
        MatrixXf r(3, 3);
        r.setZero();

        r(0,0) = q.w() * q.w() + q.x() * q.x() - q.y() * q.y() - q.z() * q.z(), r(0,1) = 2 * (q.x() * q.y() - q.w() * q.z()), r(0,2) = 2 * (q.x() * q.z() + q.w() * q.y());
        r(1,0) = 2 * (q.x() * q.y() + q.w() * q.z()), r(1,1) = q.w() * q.w() - q.x() * q.x() + q.y() * q.y() - q.z() * q.z(), r(1,2) = 2 * (q.y() * q.z() - q.w() * q.x());
        r(2,0) = 2 * (q.x() * q.z() - q.w() * q.y()), r(2,1) = 2 * (q.y() * q.z() + q.w() * q.x()), r(2,2) = q.w() * q.w() - q.x() * q.x() - q.y() * q.y() + q.z() * q.z();

        return r;
    }

    static Vector3f quaternion2rotvec(Quaternionf& q)//四元数转等效旋转矢量
    {
        Vector3f r;
        if (q.w() != 0)
        {
            float t = atan2(sqrt(q.x() * q.x() + q.y() * q.y() + q.z() * q.z()), q.w());
            r = q.vec() * ((2 * t) / sin(t));
        }
        else
            r = q.vec() * M_PI;
        return r;
    }


    static Quaternionf rotvec2quaternion(Vector3f& rotvec)//等效旋转矢量转四元数
    {
        Vector3f v;
        v = rotvec * 0.5;
        Quaternionf q(0,0,0,0);
        float t = v.norm();
        if(t==0)
        {
            q.w()=1;
            return q;
        }

        q.w() = cos(t);
        v = v * (sin(t) / t);
        q.x() = v(0);
        q.y() = v(1);
        q.z() = v(2);
        return q;
    }

    static MatrixXf get_Antisymmetric_Matrix(Vector3f& v)//向量转反对称矩阵
    {
        MatrixXf r(3, 3);
        r(0,0) = 0, r(0,1) = -v(2), r(0,2) = v(1);
        r(1,0) = v(2), r(1,1) = 0, r(1,2) = -v(0);
        r(2,0) = -v(1), r(2,1) = v(0), r(2,2) = 0;
        return r;
    }

    static MatrixXf rotvec2matrix(Vector3f& v)//等效旋转矢量转方向余弦阵
    {
        float t = v.norm();
        
        MatrixXf m(3, 3);
        MatrixXf temp(3,3);
        temp.setIdentity();

        if(t != 0)
        {

            m = get_Antisymmetric_Matrix(v);//v.cross(Eigen::Matrix3f::Identity());

            m = temp + m * (sin(t) / t) + m * m * ((1 - cos(t) )/ (t * t));

        }
        else
        {
            m.setIdentity();
        }

        return m;



    }

    static Quaternionf matrix2quaternion(MatrixXf& m)//方向余弦阵转四元数
    {
        MatrixXf t(4,1);
        t(0,0) = 1 + m.trace();
        t(1,0) = 1 + 2 * m(0,0) - m.trace();
        t(2,0) = 1 + 2 * m(1,1) - m.trace();
        t(3,0) = 1 + 2 * m(2,2) - m.trace();


        Quaternionf q;


        if (t(0,0) >= t(1,0) && t(0,0) >= t(2,0) && t(0,0) >= t(3,0))
        {
            q.w() = sqrt(t(0,0)) / 2; q.x() = (m(2,1) - m(1,2)) / (4 * q.w());
            q.y() = (m(0,2) - m(2,0)) / (4 * q.w()); q.z() = (m(1,0) - m(0,1)) / (4 * q.w());
        }
        else if (t(1,0) >= t(0,0) && t(1,0) >= t(2,0) && t(1,0) >= t(3,0))
        {
            q.x() = sqrt(t(1,0)) / 2; q.y() = (m(1,0) + m(0,1)) / (4 * q.x());
            q.z() = (m(0,2) + m(2,0)) / (4 * q.x()); q.w() = (m(2,1) - m(1,2)) / (4 * q.x());
        }
        else if (t(2,0) >= t(1,0) && t(2,0) >= t(0,0) && t(2,0) >= t(3,0))
        {
            q.y() = sqrt(t(2,0)) / 2; q.z() = (m(2,1) + m(1,2)) / (4 * q.y());
            q.w() = (m(0,2) - m(2,0)) / (4 * q.y()); q.x() = (m(0,1) + m(1,0)) / (4 * q.y());
        }
        else
        {
            q.z() = sqrt(t(3,0)) / 2; q.w() = (m(1,0) - m(0,1)) / (4 * q.z());
            q.x() = (m(0,2) + m(2,0)) / (4 * q.z()); q.y() = (m(2,1) + m(1,2)) / (4 * q.z());

        }
        return q;

    }

    static Vector3f matrix2vector3d(MatrixXf &m)//3*1列向量转向量描述
    {
        if((m.cols()==3 && m.rows()==1)||(m.rows()==3 && m.cols()==1))
        {
            Vector3f r;
            r(0)=m(0,0);
            r(1)=m(1,0);
            r(2)=m(2,0);
            return r;
        }
        Vector3f r;
        r.setZero();
        return r;
    }

    static MatrixXf vectorcrossmatrix(Vector3f &v)//求反对称矩阵
    {
        MatrixXf m(3,3);
        m.setZero();
        m(0,1) = -v(2);m(0,2) = v(1);
                        m(1,2) = -v(0);
        m(1,0) = v(2);
        m(2,0) = -v(1);m(2,1) = v(0);
        return m;
    }

    static float deg2rad(float t)//°转弧度制
    {
        return  (t / 180)*M_PI;

    }


    static float rad2deg(float t)//弧度制转角度
    {
        return t * 180 / M_PI;
    }

    static float get_g(float t)//计算重力加速度与时间间隔的乘积
    {
        return 9.7803267715 * t;
    }


};







#endif /* IMU_TRANSFORMER_H_ */
