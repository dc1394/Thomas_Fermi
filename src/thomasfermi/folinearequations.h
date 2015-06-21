/*! \file folinearequations.h
	\brief 一次要素に対して連立方程式を解くクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _FOLENEAREQUATIONS_H_
#define _FOLENEAREQUATIONS_H_

#pragma once

#include "linearequations.h"

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			一次要素に対して連立方程式を解くクラス
		*/
		class FOLinear_equations final : public Linear_equations {
			// #region コンストラクタ・デストラクタ

		public:
			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param res 連立方程式の要素
			*/
			explicit FOLinear_equations(FEM::resultmap const & res);

			//! A destructor.
			/*!
				デフォルトデストラクタ
			*/
			~FOLinear_equations() noexcept override = default;

			// #endregion コンストラクタ・デストラクタ 

			// #region publicメンバ関数

			//! A public member function.
			/*!
				連立一次方程式の解を求める
				\return 連立一次方程式の解
			*/
			FEM::dmklvector LEsolver() override;

			//! A public member function (pure virtual function).
			/*!
				初期化する
				\param b 対象のベクトルb
			*/
			void reset(FEM::dmklvector const & b) override;

			// #endregion publicメンバ関数

			// #region 禁止されたコンストラクタ・メンバ関数

			//! A private constructor (deleted).
			/*!
				デフォルトコンストラクタ（禁止）
			*/
			FOLinear_equations() = delete;

			//! A private copy constructor (deleted).
			/*!
				コピーコンストラクタ（禁止）
				\param コピー元のオブジェクト（未使用）
			*/
			FOLinear_equations(FOLinear_equations const &) = delete;

			//! A private member function (deleted).
			/*!
				operator=()の宣言（禁止）
				\param コピー元のオブジェクト（未使用）
				\return コピー元のオブジェクト
			*/
			FOLinear_equations & operator=(FOLinear_equations const &) = delete;

			// #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif	// _FOLENEAREQUATIONS_H_
