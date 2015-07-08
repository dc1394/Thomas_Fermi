/*! \file linearequations.cpp
	\brief 連立方程式を解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "fem.h"
#include "linearequations.h"
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

		Linear_equations::Linear_equations(FEM::resulttuple const & res) :
            a0_(std::get<0>(res)),
            a0back_(a0_),
            a1_(std::get<1>(res)),
            a1back_(a1_),
            a2_(std::get<2>(res)),
            b_(std::get<3>(res)),
            n_(std::get<0>(res).size())
		{
		}

		// #endregion コンストラクタ

		// #region publicメンバ関数 

		void Linear_equations::bound(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
		{
			b_[i_bc_nonzero[0] + 1] -= v_bc_nonzero[0] * a2_[i_bc_nonzero[0]];
			b_[i_bc_nonzero[1] - 1] -= v_bc_nonzero[1] * a2_[i_bc_nonzero[1] - 1];

			for (auto i = 0U; i < n_bc_nonzero; i++) {
				b_[i_bc_nonzero[i]] = v_bc_nonzero[i];
			}

			for (auto i = 0U; i < n_bc_given; i++) {
				a1_[i_bc_given[i]] = 1.0;
				a2_[i_bc_given[i] - i] = 0.0;
			}
		}
		
		Linear_equations::dvector Linear_equations::LEsolver()
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
			} else if (info < 0) {
				auto const str = (boost::format("%d-th argument has illegal value") % std::abs(info)).str();

				throw std::invalid_argument(str);
			}

			return std::move(b_);
		}

		void Linear_equations::reset(dvector const & b)
		{
			a0_ = a0back_;
			a1_ = a1back_;
			b_ = b;
		}
	}
}
