/*! \file folinearequations.cpp
	\brief 一次要素に対して連立方程式を解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "folinearequations.h"
#include <cmath>				// for std::abs
#include <utility>				// for std::move
#include <stdexcept>			// for std::logic_error
#include <boost/cast.hpp>		// for boost::numeric_cast
#include <boost/format.hpp>		// for boost::format
#include <mkl.h>				// for LAPACKE_dptsv

///* C++の複素数型を使う */
//#define LAPACK_COMPLEX_CPP
///* 関数のシンボル名は小文字 */
//#define LAPACK_NAME_PATTERN_LC
//#include "lapacke_config.h"
//#include "lapacke.h"
//
//#pragma comment(lib, "liblapacke.lib")

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ

		FOLinear_equations::FOLinear_equations(FEM::resultmap const & res)
			: Linear_equations(res)
		{
		}

		// #region コンストラクタ

		// #region publicメンバ関数

		FEM::dmklvector FOLinear_equations::LEsolver()
		{
			auto const n = boost::numeric_cast<lapack_int>(n_);
			auto const info = LAPACKE_dptsv(
				LAPACK_COL_MAJOR,
				n,
				1,
				a0_.data(),
				a1_.data(),
				b_.data(),
				n);

			if (info > 0) {
				throw std::logic_error("U is singular");
			}
			else if (info < 0) {
				auto const str = (boost::format("%d-th argument has illegal value") % std::abs(info)).str();

				throw std::invalid_argument(str);
			}

			return std::move(b_);
		}

		void FOLinear_equations::reset(FEM::dmklvector const & b)
		{
			a0_ = a0back_;
			a1_ = a1back_;
			b_ = b;
		}

		// #endregion publicメンバ関数
	}
}
