/*! \file fem.h
	\brief 有限要素法のクラスの宣言

	Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _FEM_H_
#define _FEM_H_

#pragma once

#ifdef _MSC_VER
    #pragma warning(disable : 4819)
    #define _SCL_SECURE_NO_WARNINGS
#endif

#include "beta.h"
#include "gausslegendre/gausslegendre.h"
#include "mkl_allocator.h"
#include "utility/property.h"
#include <functional>					// for std::function
#include <memory>						// for std::unique_ptr, std::shared_ptr
#include <tuple>						// for std::tuple
#include <vector>						// for std::vector
#include <boost/multi_array.hpp>		// for boost::multi_array

namespace thomasfermi {
	namespace femall {
		using namespace utility;

		//! A class.
		/*!
			有限要素法のクラス
		*/
		class FEM {
			// #region 型エイリアス

		public:
			using dvector = std::vector<double>;

			using dmklvector = std::vector<double, mkl_allocator<double>>;

		private:
			using dmatrix = boost::multi_array<double, 2>;

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
			FEM(dvector && beta, dvector const & coords, std::size_t nint, bool usesimd, bool usetbb);

			//! A destructor.
			/*!
				デストラクタ
			*/
			virtual ~FEM() = default;
			
			// #endregion コンストラクタ・デストラクタ 

			// #region publicメンバ関数

			//! A public member function (constant).
			/*!
				結果を返す関数
				\return 結果を集めたtuple
			*/
			std::tuple<dmklvector, dmklvector, dmklvector> createresult() const;
	
			//! A public member function.
			/*!
				βの状態をリセットする
				\param beta 対象のβ
			*/
			void reset(dvector const & beta);
			
			//! A public member function.
			/*!
				
			*/
			void stiff();

			//! A public member function.
			/*!

			*/
			void stiff2();

			// #endregion publicメンバ関数

			// #region protectedメンバ関数

		protected:
			//! A protected member function.
			/*!
				\param beta 新しいβ
			*/
			void initialize();

			// #endregion protectedメンバ関数

			// #region privateメンバ関数

		private:
			//! A private member function.
			/*!
				\param ielem
			*/
			void amerge(std::size_t ielem);

			//! A private member function.
			/*!
				\param ielem
			*/
			void element(std::size_t ielem);

			//! A private member function (pure virtual function).
			/*!
				\param ielem
			*/
			virtual dvector getc(std::size_t ielem) const = 0;
			
			//! A private member function (pure virtual function).
			/*!

			*/
			virtual dvector getdndr() const = 0;
						
			// #endregion 

			// #region プロパティ

		public:
			//! A property.
			/*!

			*/
			Property<FEM::dmklvector> const B;

			//! A property.
			/*!

			*/
			Property<std::size_t> const Nnode;

			//! A property.
			/*!
				
			*/
			Property<std::size_t> const Ntnoel;

			//! A property.
			/*!
				βのスマートポインタへのプロパティ
			*/
			Property<std::shared_ptr<Beta>> const PBeta;

			// #endregion プロパティ

			// #region メンバ変数
			
		protected:
			//! A protected member variable (constant).
			/*!
			*/
			std::size_t const nnode_;

		private:
			//! A private member variable.
			/*!
			*/
			dmklvector a1_;

			//! A private member variable.
			/*!
			*/
			dmklvector a2_;

			//! A private member variable.
			/*!
			*/
			dmklvector b_;

			//! A private member variable (constant).
			/*!
				関数β
			*/
			dvector const beta_;
			
			//! A private member variable.
			/*!
			*/
			std::unique_ptr<dmatrix> pastiff_;

		protected:
			//! A protected member variable (constant).
			/*!
			*/
			dvector const coords_;

			//! A protected member variable.
			/*!
				double func(double)の形の関数オブジェクト
			*/
			std::function<double(double)> func_;

			//! A protected member variable.
			/*!
				Gauss-Legendre積分のオブジェクト
			*/
			gausslegendre::Gauss_Legendre gl_;

			//! A protected member variable.
			/*!
			*/
			std::size_t nelem_;
			
			//! A protected member variable (constant).
			/*!
			*/
			std::size_t const nint_;

			//! A protected member variable.
			/*!
			*/
			std::size_t ntnoel_;
			
			//! A protected member variable.
			/*!
				βオブジェクトへのスマートポインタ
			*/
			std::shared_ptr<Beta> pbeta_;
			
			//! A protected member variable.
			/*!
			*/
			std::unique_ptr<boost::multi_array<std::size_t, 2>> plnods_;

			//! A protected member variable.
			/*!
				SIMDを使用するかどうか
			*/
			bool const usesimd_;
			
			//! A protected member variable.
			/*!
				TBBを使用するかどうか
			*/
			bool const usetbb_;

			// #region 禁止されたコンストラクタ・メンバ関数

		private:
			//! A private constructor (deleted).
			/*!
				デフォルトコンストラクタ（禁止）
			*/
			FEM() = delete;

			//! A private copy constructor (deleted).
			/*!
				コピーコンストラクタ（禁止）
				\param コピー元のオブジェクト（未使用）
			*/
			FEM(FEM const &) = delete;

			//! A private member function (deleted).
			/*!
				operator=()の宣言（禁止）
				\param コピー元のオブジェクト（未使用）
				\return コピー元のオブジェクト
			*/
			FEM & operator=(FEM const &) = delete;

			// #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif  // _FEM_H_
