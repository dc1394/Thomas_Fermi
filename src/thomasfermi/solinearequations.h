/*! \file folinearequations.h
	\brief 二次要素に対して連立方程式を解くクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _SOLENEAREQUATIONS_H_
#define _SOLENEAREQUATIONS_H_

#pragma once

#include "linearequations.h"

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			一次要素に対して連立方程式を解くクラス
		*/
		class SOLinear_equations final : public Linear_equations {
			// #region コンストラクタ・デストラクタ

		public:
			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param res 連立方程式の要素
			*/
			explicit SOLinear_equations(FEM::resultmap const & res);

			//! A destructor.
			/*!
				デフォルトデストラクタ
			*/
			~SOLinear_equations() noexcept override = default;

			// #endregion コンストラクタ・デストラクタ 

			// #region publicメンバ関数

			//! A public member function.
			/*!
				境界条件を定める
				\param n_bc_given 既知量の数
				\param i_bc_given 既知量のインデックス
				\param n_bc_nonzero 非零の既知量の数
				\param i_bc_nonzero 非零の既知量のインデックス
				\param v_bc_nonzero 非零の既知量
			*/
			void bound(std::size_t n_bc_given, Linear_equations::sivector const & i_bc_given, std::size_t n_bc_nonzero, Linear_equations::sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero) override;

			//! A public member function.
			/*!
				連立一次方程式の解を求める
				\return 連立一次方程式の解
			*/
			FEM::dmklvector LEsolver() override;

			//! A public member function (virtual function).
			/*!
				ベクトルbを初期化する
				\param b 対象のベクトルb
			*/
			void reset(FEM::dmklvector const & b) override;
			
			// #endregion publicメンバ関数

			// #region メンバ変数

		private:
			//! A private member variable.
			/*!
				ベクトルa2
			*/
			FEM::dmklvector a2_;

			//! A private member variable (constant).
			/*!
				ベクトルa2の複製
			*/
			FEM::dmklvector const a2back_;
			
			// #endregion メンバ変数

			// #region 禁止されたコンストラクタ・メンバ関数

			//! A private constructor (deleted).
			/*!
				デフォルトコンストラクタ（禁止）
			*/
			SOLinear_equations() = delete;

			//! A private copy constructor (deleted).
			/*!
				コピーコンストラクタ（禁止）
				\param コピー元のオブジェクト（未使用）
			*/
			SOLinear_equations(SOLinear_equations const &) = delete;

			//! A private member function (deleted).
			/*!
				operator=()の宣言（禁止）
				\param コピー元のオブジェクト（未使用）
				\return コピー元のオブジェクト
			*/
			SOLinear_equations & operator=(SOLinear_equations const &) = delete;

			// #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif	// _FOLENEAREQUATIONS_H_
