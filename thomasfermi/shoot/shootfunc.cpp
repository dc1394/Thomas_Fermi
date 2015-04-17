/*! \file shootfunc.h
    \brief y(x)の初期関数y0(x)の、原点に近い点xにおける関数値および微分値と、
           適合点xfにおける関数値および微分値を求めるクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "shootfunc.h"
#include <cmath>						// for std::sqrt
#include <utility>						// for std::move
#include <boost/range/algorithm.hpp>	// for boost::copy

namespace thomasfermi {
	namespace shoot {
		shootfunc::state_type shootfunc::load1(double x1, double v1)
		{
			shootfunc::state_type y;

			//y[0] = 1.0 + v1[0] * x1 + 4.0 / 3.0 * x1 * std::sqrt(x1) + 0.4 * v1[0] * x1 * x1 * std::sqrt(x1) + 1.0 / 3.0 * x1 * x1 * x1;
			y[0] = (((1.0 / 3.0 * x1 + 0.4 * v1 * std::sqrt(x1)) * x1) + 4.0 / 3.0 * std::sqrt(x1) + v1) * x1 + 1.0;
			//y[1] = v1[0] + 2.0 * std::sqrt(x1) + v1[0] * x1 * std::sqrt(x1) + x1 * x1 + 0.15 * v1[0] * x1 * x1 * std::sqrt(x1);
			y[1] = ((0.15 * v1 * std::sqrt(x1) + 1.0) * x1 + v1 * std::sqrt(x1)) * x1 + 2.0 * std::sqrt(x1) + v1;

			return std::move(y);
		}

		shootfunc::dblasvector shootfunc::score(shootfunc::state_type const & y)
		{
            shootfunc::dblasvector f(shootfunc::NVAR);
			boost::copy(y, f.begin());

			return std::move(f);
		}

		void shootfunc::rhs(const shootfunc::state_type & y, shootfunc::state_type & dydx, const double x)
		{
			dydx[0] = y[1];
			dydx[1] = y[0] * std::sqrt(y[0] / x);
		}
	}
}
