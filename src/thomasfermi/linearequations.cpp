/*! \file linearequations.cpp
\brief 連立方程式を解くクラスの実装

Copyright ©  2015 @dc1394 All Rights Reserved.
This software is released under the BSD 2-Clause License.
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

		template <>
		void Linear_equations::bound<Element::First>(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
		{
			b_[i_bc_nonzero[0] + 1] -= v_bc_nonzero[0] * a1_[i_bc_nonzero[0]];
			b_[i_bc_nonzero[1] - 1] -= v_bc_nonzero[1] * a1_[i_bc_nonzero[1] - 1];

			for (auto i = 0U; i < n_bc_nonzero; i++) {
				b_[i_bc_nonzero[i]] = v_bc_nonzero[i];
			}

			for (auto i = 0U; i < n_bc_given; i++) {
				a0_[i_bc_given[i]] = 1.0;
				a1_[i_bc_given[i] - i] = 0.0;
			}
		}

		template <>
		void Linear_equations::bound<Element::Second>(std::size_t n_bc_given, Linear_equations::sivector const & i_bc_given, std::size_t n_bc_nonzero, Linear_equations::sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
		{
			// Dirichlet boundary condition
			Linear_equations::bound<Element::First>(n_bc_given, i_bc_given, n_bc_nonzero, i_bc_nonzero, v_bc_nonzero);

			b_[i_bc_nonzero[0] + 2] -= v_bc_nonzero[0] * a2_[i_bc_nonzero[0]];
			b_[i_bc_nonzero[1] - 2] -= v_bc_nonzero[1] * a2_[i_bc_nonzero[1] - 2];

			a2_[i_bc_given[0]] = 0.0;
			a2_[i_bc_given[1] - 2] = 0.0;
		}

		template <>
		Linear_equations::dvector Linear_equations::LEsolver<Element::First>()
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

		template <>
		FEM::dmklvector Linear_equations::LEsolver<Element::Second>()
		{
			auto const n = boost::numeric_cast<lapack_int>(n_);

			// 係数行列の帯の中にある対角線より上の部分の個数
			lapack_int const kd = 2;

			// 行列{B}の列数。通常通り1
			lapack_int const nrhs = 1;

			// 配列ABの1次元目の大きさ（= KD + 1）
			auto const nb = kd + 1;

			// 係数行列の帯の外を省略して詰め込んだ2次元配列
			// ピボッティングありのLU分解を行うために(KD + 1)× N必要
			FEM::dmklvector ab(nb * n);

			for (auto i = 0; i < n; i++) {
				for (auto j = i; j <= i + 2; j++) {
					if (j == i) {
						ab[(j)* nb + (kd + i - j)] = a0_[i];
					}
					else if (j == i + 1 && j < n - 1) {
						ab[(j)* nb + (kd + i - j)] = a1_[i];
					}
					else if (j == i + 2 && j < n - 2) {
						ab[(j)* nb + (kd + i - j)] = a2_[i];
					}
				}
			}

			auto const info = LAPACKE_dpbsv(
				LAPACK_COL_MAJOR,				// 行優先か列優先か
				'U',							// 上三角要素を使う場合
				n,								// 線形方程式の数（行列Aの次数）
				kd,								// 係数行列の帯の中にある対角線より上の部分の個数
				nrhs,							// 行列{B}の列数。通常通り1
				ab.data(),						// 係数行列(input)，コレスキー分解の結果(output)
				nb,								// 配列ABの1次元目の大きさ（=KD+1） 
				b_.data(),						// 方程式の右辺(input)，方程式の解(output)
				n);								// 行列Bの1次元目の大きさ（=N）

			if (info > 0) {
				throw std::logic_error("U is singular");
			}
			else if (info < 0) {
				auto const str = (boost::format("%d-th argument has illegal value") % std::abs(info)).str();

				throw std::invalid_argument(str);
			}

			return std::move(b_);
		}

		void Linear_equations::reset(Linear_equations::dvector const & b)
		{
			a0_ = a0back_;
			a1_ = a1back_;
			b_ = b;
		}

		// #endregion publicメンバ関数

		// #region templateメンバ関数の実体化

		template <>
		void Linear_equations::bound<Element::First>(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero);

		template <>
		void Linear_equations::bound<Element::Second>(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero);

		template <>
		Linear_equations::dvector Linear_equations::LEsolver<Element::First>();

		template <>
		Linear_equations::dvector Linear_equations::LEsolver<Element::Second>();

		// #endregion templateメンバ関数の実体化
	}
}
