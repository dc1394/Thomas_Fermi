/*! \file linearequations.cpp
	\brief 連立方程式を解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "linearequations.h"

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ

		Linear_equations::Linear_equations(FEM::resultmap const & res) :
			a0_(res.at("a0")),
			a0back_(a0_),
			a1_(res.at("a1")),
			a1back_(a1_),
			b_(res.at("b")),
			n_(a0_.size())
		{
		}

		// #endregion コンストラクタ

		// #region publicメンバ関数
		
		void Linear_equations::bound(std::size_t n_bc_given, Linear_equations::sivector const & i_bc_given, std::size_t n_bc_nonzero, Linear_equations::sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
		{
			// Dirichlet boundary condition
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

		// #endregion publicメンバ関数 
	}
}
