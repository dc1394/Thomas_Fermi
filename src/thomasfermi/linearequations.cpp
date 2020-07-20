/*! \file linearequations.cpp
    \brief 連立方程式を解くクラスの実装
    Copyright © 2015-2019 @dc1394 All Rights Reserved.
    
    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "fem.h"
#include "linearequations.h"

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        Linear_equations::Linear_equations(FEM::resulttuple const & res)
            :   a0_(std::get<0>(res)),
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

        void Linear_equations::reset(std::vector<double> const & b)
        {
            a0_ = a0back_;
            a1_ = a1back_;
            b_ = b;
        }

        // #endregion publicメンバ関数
    }
}
