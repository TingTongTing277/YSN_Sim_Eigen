# Sim_Eigen

## ��Ŀ����
����һ���򵥵ľ����㷨�⣬�����˾�����ά��������Ԫ����ʹ�ñ�׼C++��д��
��Ҫ����һЩ�ϲ��MCU����������������֤�㷨��
## ���󻷾�
��׼C++������֧��iostream��stdexcept
##ʹ�÷���
�������ʹ��Eigenһ��ʹ���������磺
```cpp
#include <iostream>
#include "Math_Header.h"
MatrixXf m(3,3);
m << 1,2,3,
    4,5,6,
    7,8,9;
std::cout << m << std::endl;
```

## Ŀ¼�ṹ����
������ ReadMe.md           // �����ĵ�

������ Block.h              	// �ֿ������

������ Matrix.h           	// ������ʵ��

������ Vector.h           	// ������ʵ��

������ Quaterniond.h     	// ��Ԫ����ʵ��

������ Transformer.h     	// ��̬�໥�任

������ Math_Header.h     // ��ͷ�ļ�

## �汾���ݸ���
V1.0.0: ��һ�θ��� 2025.5.12
