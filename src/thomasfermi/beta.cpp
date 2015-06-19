﻿/*! \file beta.cpp
	\brief β(x)を計算するクラスの実装

	Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "beta.h"

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ

		Beta::Beta(std::vector<double> const & xvec, std::vector<double> const & yvec) :
			size_(xvec.size()),
			xvec_(xvec),
			yvec_(yvec)
		{
		}
		
		// #endregion コンストラクタ
	}
}
