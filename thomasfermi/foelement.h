/*! \file foelement.h
	\brief 一次要素のクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _FOELEMENT_H_
#define _FOELEMENT_H_

#pragma once

#include "fem.h"

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			一次要素のクラス
		*/
		class FOElement final : public FEM {
			// #region コンストラクタ・デストラクタ

		public:
			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param beta
				\param coords
				\param nint
				\param usesimd SIMDを利用するかどうか
				\param usetbb TBBを使用するかどうか
			*/
			FOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usesimd, bool usetbb);

			//! A destructor.
			/*!
				デストラクタ
			*/
			virtual ~FOElement() = default;

			// #endregion コンストラクタ・デストラクタ 

			// #region メンバ関数

		private:
			//! A private member function (constant).
			/*!
				dn/drを返す関数
				\return dn/dr
			*/
			dvector getdndr() const;

			//! A private member function (constant).
			/*!
				cを返す関数
				\param ielem 要素
				\return c
			*/
			dvector getc(std::size_t ielem) const;

			// #endregion メンバ関数

			// #region メンバ変数
			
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
			std::function<double (double)> const N1_;
			
			//! A private member variable (constant).
			/*!
				関数オブジェクト
			*/
			std::function<double (double)> const N2_;

			// #region 禁止されたコンストラクタ・メンバ関数

			//! A private constructor (deleted).
			/*!
				デフォルトコンストラクタ（禁止）
			*/
			FOElement() = delete;

			//! A private copy constructor (deleted).
			/*!
				コピーコンストラクタ（禁止）
				\param コピー元のオブジェクト（未使用）
			*/
			FOElement(FOElement const &) = delete;
			
			//! A private member function (deleted).
			/*!
				operator=()の宣言（禁止）
				\param コピー元のオブジェクト（未使用）
				\return コピー元のオブジェクト
			*/
			FOElement & operator=(FOElement const &) = delete;
			
			// #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif	// _FOELEMENT_H_
