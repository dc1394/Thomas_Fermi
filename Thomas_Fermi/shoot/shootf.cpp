#include "shootf.h"
#include <iterator>
#include <cmath>
#include <cstdint>
#include <boost/cast.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/foreach.hpp>
#endif

namespace Thomas_Fermi {
	namespace shoot {
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
		const double shootf::EPS = 1.0E-14;
#endif
		shootf::result_type shootf::operator()(double x1, double x2, double xf)
		{
			BOOST_ASSERT(x1 < xf);
			BOOST_ASSERT(x2 > xf);

			using namespace boost::numeric::odeint;

			bulirsch_stoer<shootfunc::state_type> stepper(eps_, eps_);

			shootfunc::state_type y1 = load1_(x1, v1_);					// 最良の仮の値v1_でx1からxfまで解いていく
			integrate_const(stepper, shootfunc::rhs, y1, x1, xf, dx_);
			const shootfunc::dblasvector f1(score_(y1));
		
			shootfunc::state_type y2 = load2_(x2, v2_);					// 最良の仮の値v2_でx2からxfまで解いていく			
			integrate_const(stepper, shootfunc::rhs, y2, x2, xf, - dx_);
			const shootfunc::dblasvector f2(score_(y2));

			boost::numeric::ublas::matrix<double> dfdv(shootfunc::NVAR, shootfunc::NVAR);
			std::int32_t j = 0;
			for (std::size_t iv = 0; iv < shootfunc::N2; iv++, j++) {	// x1で用いる境界条件を変える
				const double sav = v1_[iv];
				v1_[iv] += delv1_[iv];

				shootfunc::state_type y = load1_(x1, v1_);
				integrate_const(stepper, shootfunc::rhs, y, x1, xf, dx_);
				const shootfunc::dblasvector f = score_(y);

				for (std::size_t i = 0; i < shootfunc::NVAR; i++)		// NVAR個の合致条件にある偏微分を数値的に計算
					dfdv(i, j) = (f[i] - f1[i]) / delv1_[iv];
			

				v1_[iv] = sav;											// 境界におけるパラメータを格納
			}

			for (std::size_t iv = 0; iv < shootfunc::N2; iv++, j++) {	// 次にx2で用いる境界条件を変える		
				const double sav = v2_[iv];
				v2_[iv] += delv2_[iv];

				shootfunc::state_type y = load2_(x2, v2_);
				integrate_const(stepper, shootfunc::rhs, y, x2, xf, - dx_);
				const shootfunc::dblasvector f = score_(y);

				for (std::size_t i = 0; i < shootfunc::NVAR; i++)
					dfdv(i, j) = (f2[i] - f[i]) / delv2_[iv];
			
				v2_[iv] = sav;
			}

			shootfunc::dblasvector f(shootfunc::NVAR), ff(shootfunc::NVAR);
			for (std::size_t i = 0; i < shootfunc::NVAR; i++) {
				f[i] = f1[i] - f2[i];
				ff[i] = - f[i];
			}

			boost::numeric::ublas::permutation_matrix<> pm(dfdv.size1());
			boost::numeric::ublas::lu_factorize(dfdv, pm);				// 自由パラメータに対する増分を求める
			boost::numeric::ublas::lu_substitute(dfdv, pm, ff);

			j = 0;
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			for (double & v : v1_)										// x1の境界でのパラメータ値の増分
				v += ff[j];

			for (double & v : v2_)										// x2の境界でのパラメータ値の増分
				v += ff[++j];
#else
			BOOST_FOREACH (double & v, v1_)								// x1の境界でのパラメータ値の増分
				v += ff[j];
			
			BOOST_FOREACH (double & v, v2_)								// x2の境界でのパラメータ値の増分
				v += ff[++j];
#endif

			y1 = load1_(x1, v1_);										
			dvector res1;
			res1.reserve(boost::numeric_cast<std::size_t>((xf - x1) / dx_) + 2);
			// 得られた条件でx1...dxまで微分方程式を解く
			integrate_const(stepper, shootfunc::rhs, y1, x1, dx_, dx_ - x1, [&res1](const shootfunc::state_type & y, const double x)
			{ res1.push_back(y[0]); });									// x1...dxの結果を得る
			res1.pop_back();

			// 得られた条件でdx...xfまで微分方程式を解く
			integrate_const(stepper, shootfunc::rhs, y1, dx_, xf + x1, dx_, [&res1](const shootfunc::state_type & y, const double x)
			{ res1.push_back(y[0]); });									// dx...xf + x1の結果を得る
			
			// 得られた条件でx2...xfまで微分方程式を解く
			y2 = load2_(x2, v2_);								
			dvector res2;
			res2.reserve(boost::numeric_cast<std::size_t>((x2 - xf) / dx_) + 1);
			integrate_const(stepper, shootfunc::rhs, y2, x2, xf - x1, - dx_, [&res2](const shootfunc::state_type & y, const double x)
			{ res2.push_back(y[0]); });									// x2...xf - x1の結果を得る

			return createResult(x1, xf, res1, res2);
		}

		shootf::result_type shootf::createResult(double x1, double xf,
												 const dvector & res1, const dvector & res2) const
		{
			const std::size_t size = res1.size() + res2.size() - 1;
			dvector xp, xptmp, yp;
			xp.reserve(size);
			xptmp.reserve(size - 1);
			yp.reserve(size);

			const std::size_t xfindex = boost::numeric_cast<std::size_t>(xf / dx_); 
			for (std::size_t i = 0; i < xfindex; i++)
				xptmp.push_back(i ? static_cast<double>(i) * dx_ : x1);
			for (std::size_t i = xfindex + 1; i < size; i++)
				xptmp.push_back(static_cast<double>(i) * dx_);

			yp.assign(res1.begin(), res1.end() - 1);
			const std::size_t s = res2.size();
			for (std::size_t i = 1; i < s; i++)
				yp.push_back(res2[s - i - 1]);

			BOOST_ASSERT(xptmp.size() == yp.size());

			const Spline splint(xptmp, yp);
			
			for (std::size_t i = 0; i < size; i++)
				xp.push_back(i ? static_cast<double>(i) * dx_ : x1);

			auto iter2 = yp.begin();
			std::advance(iter2, std::distance(iter2, boost::find_if(yp,
				[&, s](double x) { return std::fabs(x - res2[s - 2]) < EPS; })));
			yp.insert(iter2, splint(xf));

			BOOST_ASSERT(xp.size() == yp.size());

			return std::make_tuple(std::move(xp), std::move(yp), v1_);
		}
	}
}
