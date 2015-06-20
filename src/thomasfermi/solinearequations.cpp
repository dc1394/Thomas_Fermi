/*! \file solinearequations.cpp
	\brief 二次要素に対して連立方程式を解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "solinearequations.h"
#include <cmath>				// for std::abs
#include <utility>				// for std::move
#include <stdexcept>			// for std::logic_error
#include <vector>				// for std::vector
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

		SOLinear_equations::SOLinear_equations(FEM::resultmap const & res) :
			Linear_equations(res),
			a2_(res.at("a2")),
			a2back_(a2_)
		{
		}

		// #region コンストラクタ

		// #region publicメンバ関数

		void SOLinear_equations::bound(std::size_t n_bc_given, Linear_equations::sivector const & i_bc_given, std::size_t n_bc_nonzero, Linear_equations::sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
		{
			// Dirichlet boundary condition
			b_[i_bc_nonzero[0] + 1] -= v_bc_nonzero[0] * (a1_[i_bc_nonzero[0]] + a2_[i_bc_nonzero[0]]);
			b_[i_bc_nonzero[1] - 1] -= v_bc_nonzero[1] * (a1_[i_bc_nonzero[1] - 1] + a2_[i_bc_nonzero[1] - 1]);

			for (auto i = 0U; i < n_bc_nonzero; i++) {
				b_[i_bc_nonzero[i]] = v_bc_nonzero[i];
			}

			for (auto i = 0U; i < n_bc_given; i++) {
				a0_[i_bc_given[i]] = 1.0;
				a1_[i_bc_given[i] - i] = 0.0;
			}

			a2_[i_bc_given[0]] = 0.0;
			a2_[i_bc_given[1] - 2] = 0.0;
		}

		FEM::dmklvector SOLinear_equations::LEsolver()
		{
			auto const n = boost::numeric_cast<lapack_int>(n_);
			
			// 係数行列の帯の中にある subdiagonals (対角線より下の部分) の個数
			lapack_int const kl = 2;

			// 係数行列の帯の中にある superdiagonals (対角線より上の部分) の個数
			lapack_int const ku = 2;

			// 行列{B}の列数。通常通り1
			lapack_int const nrhs = 1;

			// 係数行列の帯の外を省略して詰め込んだ2次元配列
			// ピボッティングありのLU分解を行うために(2KL + KU + 1)× N必要
			std::vector< double *, mkl_allocator < double * > > AB(2 * kl + ku + 1);
			for (auto & elem : AB) {
				elem = mymklalloc<double>(n);
				/*elem = reinterpret_cast<double *>(mkl_malloc(n * sizeof(double), 64));

				if (!elem) {
					throw std::bad_alloc();
				}*/
			}

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

		void SOLinear_equations::reset(FEM::dmklvector const & b)
		{
			Linear_equations::reset(b);
			a2_ = a2back_;
		}

		// #endregion publicメンバ関数
	}
}
