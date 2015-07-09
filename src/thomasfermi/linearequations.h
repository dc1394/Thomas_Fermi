/*! \file linearequations.h
	\brief 連立方程式を解くクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _LENEAREQUATIONS_H_
#define _LENEAREQUATIONS_H_

#pragma once

#ifdef _MSC_VER
#pragma warning(disable : 4819)
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "element.h"
#include "mkl_allocator.h"
#include <tuple>			// for std::tuple
#include <vector>			// for std::vexctor

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			連立方程式を解くクラス
		*/
		class Linear_equations final {
			// #region 型エイリアス

			using dvector = std::vector < double, mkl_allocator<double> > ;

			using sivector = std::vector<std::size_t>;

			// #endregion 型エイリアス

			// #region コンストラクタ・デストラクタ

		public:
			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param res 連立方程式のベクトル
			*/
			explicit Linear_equations(FEM::resulttuple const & res);

			//! A destructor.
			/*!
				デストラクタ
			*/
			~Linear_equations() = default;

			// #endregion コンストラクタ・デストラクタ 

			// #region publicメンバ関数

			template <Element E>
			//! A public member function.
			/*!
				境界条件を定める
				\param n_bc_given 既知量の数
				\param i_bc_given 既知量のインデックス
				\param n_bc_nonzero 非零の既知量の数
				\param i_bc_nonzero 非零の既知量のインデックス
				\param v_bc_nonzero 非零の既知量
			*/
			void bound(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero); 
			
			//! A public member function.
			/*!
				ベクトルbを初期化する
				\param b 対象のベクトルb
			*/
			void reset(dvector const & b);
			
			template <Element E>
			//! A public member function.
			/*!
				連立一次方程式の解を求める
				\return 連立一次方程式の解
			*/
			Linear_equations::dvector LEsolver();

			// #endregion publicメンバ関数

			// #region メンバ変数

		private:
            //! A private member variable.
            /*!
                ベクトルa0
            */
            dvector a0_;

            //! A private member variable (constant).
            /*!
                ベクトルa0の複製
            */
            dvector const a0back_;

            //! A private member variable.
            /*!
                ベクトルa1
            */
            dvector a1_;

            //! A private member variable (constant).
            /*!
                ベクトルa1の複製
            */
            dvector const a1back_;

            //! A private member variable.
            /*!
                ベクトルa2
            */
            dvector a2_;

			//! A private member variable.
			/*!
				ベクトルb
			*/
			dvector b_;

			//! A private member variable (constant).
			/*!
				ベクトルの要素数
			*/
			std::size_t const n_;

			// #endregion メンバ変数		
			
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
