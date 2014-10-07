#include "shootfunc.h"
#include <cmath>
#include <utility>
#include <boost/range/algorithm.hpp>

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
