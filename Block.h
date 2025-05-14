#pragma once
#ifndef MATH_BLOCK_H_
#define MATH_BLOCK_H_

#include "Matrix.h"
#include <vector>

namespace Sim_Eigen {
	template<typename Scalar>
	class Block {
	private:
		int rows_;
		int cols_;
		std::vector<Scalar*> block_;  // 用 vector 管理指针数组
		Matrix<Scalar>* parent; //指向原始矩阵
	public:
		Block(Matrix<Scalar>& m,int startRow,int startCol, int Rows, int Cols) {
			if(startRow<0 || startCol<0 || Rows<0 || Cols<0 || startRow+Rows>m.rows() || startCol+Cols>m.cols())
				throw std::runtime_error("Matrix Block over limit");
			rows_ = Rows;
			cols_ = Cols;
			parent = &m;
			for (int i = 0; i < rows_; i++) {
				block_.push_back(&(m(startRow + i, startCol)));
			}
		}
		~Block() = default;

	public:
		int rows() const { return rows_; }
		int cols() const { return cols_; }
		Scalar** data() { return block_; }
		Matrix<Scalar>* getParent() const { return parent; }//获取原始矩阵位置
		const Scalar** data() const { return block_; }

		Scalar& operator()(int rows, int cols) {
			if (rows >= rows_ || cols >= cols_ || rows < 0 || cols < 0) {
				throw std::invalid_argument("Block element access is illegal");
			}
			return block_[rows] [cols];
		}

		const Scalar& operator()(int rows, int cols) const {
			if (rows >= rows_ || cols >= cols_ || rows < 0 || cols < 0) {
				throw std::invalid_argument("Block element access is illegal");
			}
			return block_[rows] [cols];
		}

		//分块矩阵m放回母矩阵原位置
		Block<Scalar>& operator = (const Matrix<Scalar>& m) {
			if (rows_ != m.rows() || cols_ != m.cols())
				throw std::invalid_argument("Matrix and Block size is not same!");
			for (int i = 0; i < rows_; i++)
			{
				std::copy(&(m(i,0)), &(m(i,0)) + cols_, block_[i]);
			}
			return *this;
		}

		//Block转化为矩阵输出
		Matrix<Scalar> toMatrix(void) {
			Matrix<Scalar> m(rows_, cols_);
			m = *this;
			return m;
		}
	};
}
#endif
