# Sim_Eigen

## 项目介绍
这是一个简单的矩阵算法库，包含了矩阵，三维向量，四元数，使用标准C++编写，
主要用于一些较差的MCU编译器环境快速验证算法。
## 需求环境
标准C++环境，支持iostream和stdexcept
##使用方法
你可以像使用Eigen一样使用它，比如：
```cpp
#include <iostream>
#include "Math_Header.h"
MatrixXf m(3,3);
m << 1,2,3,
    4,5,6,
    7,8,9;
std::cout << m << std::endl;
```

## 目录结构描述
├── ReadMe.md           // 帮助文档

├── Block.h              	// 分块矩阵类

├── Matrix.h           	// 矩阵类实现

├── Vector.h           	// 向量类实现

├── Quaterniond.h     	// 四元数类实现

├── Transformer.h     	// 姿态相互变换

├── Math_Header.h     // 库头文件

## 版本内容更新
V1.0.0: 第一次更新 2025.5.12
