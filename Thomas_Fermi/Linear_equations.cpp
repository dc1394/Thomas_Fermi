#include "Linear_equations.h"
#include <sstream>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <boost/cast.hpp>
#include <mkl.h>
///* C++の複素数型を使う */
//#define LAPACK_COMPLEX_CPP
///* 関数のシンボル名は小文字 */
//#define LAPACK_NAME_PATTERN_LC
//#include "lapacke_config.h"
//#include "lapacke.h"
//
//#pragma comment(lib, "liblapacke.lib")

namespace Thomas_Fermi {
	namespace FEM_ALL {
		void Linear_equations::reset(const dvector & b)
		{
			a1_ = a1back_;
			a2_ = a2back_;
			b_ = b;
		}

		void Linear_equations::bound(std::size_t n_bc_given,					// 既知量の数
									 const sivector & i_bc_given,				// 既知量のインデックス
									 std::size_t n_bc_nonzero,					// 非零の既知量の数
									 const sivector & i_bc_nonzero,				// 非零の既知量のインデックス
									 const std::vector<double> & v_bc_nonzero)	// 非零の既知量
		{
			b_[i_bc_nonzero[0] + 1] -= v_bc_nonzero[0] * a2_[i_bc_nonzero[0]];
			b_[i_bc_nonzero[1] - 1] -= v_bc_nonzero[1] * a2_[i_bc_nonzero[1] - 1];

			for (std::size_t i = 0; i < n_bc_nonzero; i++)
				b_[i_bc_nonzero[i]] = v_bc_nonzero[i];

			for (std::size_t i = 0; i < n_bc_given; i++) {
				a1_[i_bc_given[i]] = 1.0;
				a2_[i_bc_given[i] - i] = 0.0;
			}
		}

		Linear_equations::dvector Linear_equations::LEsolver()
		{
			const lapack_int n = boost::numeric_cast<lapack_int>(n_);
			const std::int32_t info = LAPACKE_dptsv(LAPACK_COL_MAJOR,
													n, 1, a1_.data(), a2_.data(),
													b_.data(), n);

			if (info > 0) {
				throw std::logic_error("U is singular");
			} else if (info < 0) {
				std::ostringstream oss;
				oss << std::abs(info) << "-th argument has illegal value";

				throw std::invalid_argument(oss.str());
			}

			return std::move(b_);
		}
	}
}
