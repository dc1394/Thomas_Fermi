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

		void Linear_equations::reset(FEM::dmklvector const & b)
		{
			a0_ = a0back_;
			a1_ = a1back_;
			b_ = b;
		}
	}
}
