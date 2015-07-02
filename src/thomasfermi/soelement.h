/*! \file soelement.h
	\brief 二次要素のクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _SOELEMENT_H_
#define _SOELEMENT_H_

#pragma once

#include "fem.h"

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			二次要素のクラス
		*/
		class SOElement final : public FEM {
			// #region コンストラクタ・デストラクタ

		public:
			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param beta
				\param coords
				\param nint
				\param usecilk Cilk Plusを使用するかどうか
			*/
			SOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk);

			//! A destructor.
			/*!
				デストラクタ
			*/
			~SOElement() noexcept override = default;

			// #endregion コンストラクタ・デストラクタ 

			// #region publicメンバ関数

			//! A public member function (constant).
			/*!
				結果を返す関数
				\return 結果を集めたboost::container::flat_map
			*/
			std::tuple<FEM::dmklvector, FEM::dmklvector, FEM::dmklvector, FEM::dmklvector> createresult() const override;

			// #endregion publicメンバ関数

			// #region privateメンバ関数

		private:
			//! A private member function.
			/*!
				a0_、a1_とa2_にastiff_を足し込む関数
				\param ielem
			*/
			void amerge(std::size_t ielem) override;

			//! A private member function.
			/*!
				小行列の要素を生成する
				\param ielem
			*/
			void element(std::size_t ielem) override;

			//! A private member function (constant).
			/*!
				dn/drを返す関数
				\param r rの値
				\return dn/dr
			*/
			dvector getdndr(double r) const override;

			//! A private member function (constant).
			/*!
				cを返す関数
				\param ielem 要素
				\return c
			*/
			dvector getc(std::size_t ielem) const override;

			// #endregion メンバ関数

			// #region メンバ変数

			//! A private member variable.
			/*!
				連立方程式Ax = Bの行列Aの三番目の要素
			*/
			dmklvector a2_;

			//! A private member variable (constant).
			/*!
				関数オブジェクト
			*/
			std::function<double(double, double, std::size_t)> const fun1_;

			//! A private member variable (constant).
			/*!
				関数オブジェクト
			*/
			std::function<double(double, double, std::size_t)> const fun2_;

			//! A private member variable (constant).
			/*!
				関数オブジェクト
			*/
			std::function<double(double, double, std::size_t)> const fun3_;

			//! A private member variable.
			/*!
				形状関数その1を格納する関数オブジェクト
			*/
			std::function<double(double)> N1_;

			//! A private member variable.
			/*!
				形状関数その2を格納する関数オブジェクト
			*/
			std::function<double(double)> N2_;

			//! A private member variable.
			/*!
				形状関数その3を格納する関数オブジェクト
			*/
			std::function<double(double)> N3_;
			
			// #region 禁止されたコンストラクタ・メンバ関数

			//! A private constructor (deleted).
			/*!
				デフォルトコンストラクタ（禁止）
			*/
			SOElement() = delete;

			//! A private copy constructor (deleted).
			/*!
				コピーコンストラクタ（禁止）
				\param コピー元のオブジェクト（未使用）
			*/
			SOElement(SOElement const &) = delete;

			//! A private member function (deleted).
			/*!
				operator=()の宣言（禁止）
				\param コピー元のオブジェクト（未使用）
				\return コピー元のオブジェクト
			*/
			SOElement & operator=(SOElement const &) = delete;

			// #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif	// _SOELEMENT_H_
