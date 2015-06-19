/*! \file linearequations.h
	\brief 連立方程式を解くクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _LENEAREQUATIONS_H_
#define _LENEAREQUATIONS_H_

#pragma once

#include "fem.h"
#include "mkl_allocator.h"
#include <vector>			// for std::vexctor

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			連立方程式を解くクラス
		*/
		class Linear_equations {
			// #region 型エイリアス

		protected:
			using sivector = std::vector<std::size_t>;

			// #endregion 型エイリアス

			// #region コンストラクタ・デストラクタ

		public:
			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param res 連立方程式の要素
			*/
			explicit Linear_equations(FEM::resultmap const & res);

			//! A destructor.
			/*!
				デフォルトデストラクタ
			*/
			virtual ~Linear_equations() noexcept = default;

			// #endregion コンストラクタ・デストラクタ 

			// #region publicメンバ関数

			//! A public member function (pure virtual function).
			/*!
				境界条件を定める
				\param n_bc_given 既知量の数
				\param i_bc_given 既知量のインデックス
				\param n_bc_nonzero 非零の既知量の数
				\param i_bc_nonzero 非零の既知量のインデックス
				\param v_bc_nonzero 非零の既知量
			*/
			virtual void bound(std::size_t n_bc_given, Linear_equations::sivector const & i_bc_given, std::size_t n_bc_nonzero, Linear_equations::sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero) = 0;
			
			//! A public member function (pure virtual function).
			/*!
				連立一次方程式の解を求める
				\return 連立一次方程式の解
			*/
			virtual FEM::dmklvector LEsolver() = 0;
			
			//! A public member function (virtual function).
			/*!
				ベクトルbを初期化する
				\param b 対象のベクトルb
			*/
			virtual void reset(FEM::dmklvector const & b);

			// #endregion publicメンバ関数

			// #region protectedメンバ変数

		protected:
			//! A private member variable.
			/*!
				ベクトルa1
			*/
			FEM::dmklvector a0_;
			
			//! A private member variable (constant).
			/*!
				ベクトルa1の複製
			*/
			FEM::dmklvector const a0back_;

			//! A private member variable.
			/*!
				ベクトルa2
			*/
			FEM::dmklvector a1_;

			//! A private member variable (constant).
			/*!
				ベクトルs2の複製
			*/
			FEM::dmklvector const a1back_;

			//! A private member variable.
			/*!
				ベクトルb
			*/
			FEM::dmklvector b_;

			//! A private member variable (constant).
			/*!
				ベクトルの要素数
			*/
			std::size_t const n_;

			// #endregion protectedメンバ変数		
			
			// #region 禁止されたコンストラクタ・メンバ関数

			//! A private constructor (deleted).
			/*!
				デフォルトコンストラクタ（禁止）
			*/
			Linear_equations() = delete;

			//! A private copy constructor (deleted).
			/*!
				コピーコンストラクタ（禁止）
				\param コピー元のオブジェクト（未使用）
			*/
			Linear_equations(Linear_equations const &) = delete;

			//! A private member function (deleted).
			/*!
				operator=()の宣言（禁止）
				\param コピー元のオブジェクト（未使用）
				\return コピー元のオブジェクト
			*/
			Linear_equations & operator=(Linear_equations const &) = delete;

			// #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif	// _LENEAREQUATIONS_H_
