/*! \file shootfunc.h
    \brief y(x)の初期関数y0(x)の、原点に近い点xにおける関数値および微分値と、
           適合点xfにおける関数値および微分値を求めるクラスの実装
    Copyright © 2014-2019 @dc1394 All Rights Reserved.

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

#include "shootfunc.h"
#include <cmath>						// for std::sqrt

namespace thomasfermi {
	namespace shoot {
		shootfunc::state_type shootfunc::load1(double x1, double v1)
		{
			shootfunc::state_type y;

			//y[0] = 1.0 + v1[0] * x1 + 4.0 / 3.0 * x1 * std::sqrt(x1) + 0.4 * v1[0] * x1 * x1 * std::sqrt(x1) + 1.0 / 3.0 * x1 * x1 * x1;
			y[0] = (((1.0 / 3.0 * x1 + 0.4 * v1 * std::sqrt(x1)) * x1) + 4.0 / 3.0 * std::sqrt(x1) + v1) * x1 + 1.0;
			//y[1] = v1[0] + 2.0 * std::sqrt(x1) + v1[0] * x1 * std::sqrt(x1) + x1 * x1 + 0.15 * v1[0] * x1 * x1 * std::sqrt(x1);
			y[1] = ((0.15 * v1 * std::sqrt(x1) + 1.0) * x1 + v1 * std::sqrt(x1)) * x1 + 2.0 * std::sqrt(x1) + v1;

			return y;
		}

		Eigen::VectorXd shootfunc::score(shootfunc::state_type const & y)
		{
			return Eigen::VectorXd::Map(y.data(), y.size());
		}

		void shootfunc::rhs(shootfunc::state_type const & y, shootfunc::state_type & dydx, double const x)
		{
			dydx[0] = y[1];
			dydx[1] = y[0] * std::sqrt(y[0] / x);
		}
	}
}
