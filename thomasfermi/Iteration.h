/*! \file iteration.h
	\brief 微分方程式を反復法で解くクラスの宣言

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _ITERATION_H_
#define _ITERATION_H_

#pragma once


#include "foelement.h"
#include "Linear_equations.h"
#include "shoot/shootfunc.h"
#include <boost/optional.hpp>	// for boost::optional;

namespace thomasfermi {
	namespace fem_all {
		class Iteration final {
			// #region 型エイリアス

		public:
			using result_type = std::tuple<FEM::dvector, std::shared_ptr<const Beta>, double>;

			// #endregion 型エイリアス

			// #region コンストラクタ・デストラクタ

			//! A constructor.
			/*!
				\param alpha 一次混合の値α
				\param dx 要素の間隔
				\param n Gauss-Legendreの積分点
				\param tol 許容誤差
				\param usesimd SIMDを使用するかどうか
				\param usetbb TBBを使用するかどうか
				\param x1 原点に近い方のxの値
				\param x2 無限遠に近い方のxの値
				\param xf 適合点
			*/
			Iteration(double alpha, double dx, std::size_t n, double eps, bool usesimd, bool usetbb, double x1, double x2, double xf);

			//! A destructor.
			/*!
				デストラクタ
			*/
			~Iteration();

			// #region コンストラクタ・デストラクタ
			
			// #region publicメンバ関数

			//! A public member function.
			/*!
				反復関数
			*/
			void Iterationloop();

			//! A public member function.
			/*!
				結果を返す関数
				\return 結果
			*/
			result_type makeresult();
			/*{
				return std::make_tuple(std::move(x_), std::move(pbeta_), v1_);
			}*/

			// #endregion publicメンバ関数

			// #region privateメンバ関数

			FEM::dvector make_beta() const;
			
			void ymix();
			
			double IterationError() const;

			// #endregion privateメンバ関数


#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
			static constexpr std::size_t N_BC_GIVEN = 2;
			static constexpr double IterationTHRESHOLD = 1.0;
			static constexpr double IterationREDUCTION = 0.15;
			const bool useSSEorAVX_;
			const bool usecilk_;
			const bool avxSupported_;
			const double TOL_;
			double alpha_;
			FEM::dvector x_;
			FEM::dmklvector y_;
			FEM::dmklvector ybefore_;
			std::unique_ptr<FEM> pfem_;
			std::shared_ptr<const Beta> pbeta_;
			boost::optional<Linear_equations> ple_;
			std::vector<std::size_t> i_bc_given_;
			FEM::dvector v_bc_nonzero_;

			double y1_;
			double y2_;
			double v1_;

			
		public:
			
				// 、無限遠に近い方、適合点、要素の間隔
			;
			~Iteration() { ple_ = boost::none; }


#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			Iteration(const Iteration &) = delete;
			Iteration & operator=(const Iteration &) = delete;
			Iteration() = delete;
#endif

		};
	}

	template <typename T>
	inline T sqr(T x)
	{ return x * x; }
}

#endif	// _ITERATION_H_
