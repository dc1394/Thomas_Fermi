/*! \file shootf.cpp
    \brief 狙い撃ち法により、y(x)を求めるクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "shootf.h"
#include <cmath>						// for std::fabs
#include <iterator>						// for std::advance, std::distance
#include <boost/assert.hpp>				// for BOOST_ASSERT
#include <boost/cast.hpp>				// for boost::numeric_cast
#include <boost/numeric/odeint.hpp>		// for boost::numeric::odeint
#include <boost/numeric/ublas/lu.hpp>	// for boost::numeric::ublas::lu_factorize, boost::numeric::ublas::lu_substitute, boost::numeric::ublas::permutation_matrix
#include <boost/range/algorithm.hpp>	// for boost::find_if			
#include <cilk/cilk.h>					// for cilk_spawn, cilk_sync

namespace thomasfermi {
	namespace shoot {
		shootf::shootf(double delv1, double delv2, double dx, double eps, loadfunctype const & load1, load2 const & l2,	scorefunctype const & score, double v1,	double v2) :
			delv1_(delv1),
			delv2_(delv2),
			dx_(dx),
			eps_(eps),
			load1_(load1),
			load2_([&l2](double v2, double x2) { return l2(v2, x2); }),
			score_(score),
			v1_(v1),
			v2_(v2)
		{
		}

		shootf::result_type shootf::operator()(double x1, double x2, double xf)
		{
			BOOST_ASSERT(x1 < xf);
			BOOST_ASSERT(x2 > xf);

			using namespace boost::numeric::odeint;

			using stepper_type = bulirsch_stoer<shootfunc::state_type>;

			// 最良の仮の値v1_でx1からxfまで解いていく
			auto y1 = load1_(x1, v1_);
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y1, x1, xf, dx_);
			auto const f1(score_(y1));

			// 最良の仮の値v2_でx2からxfまで解いていく			
			auto y2 = load2_(x2, v2_);
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y2, x2, xf, - dx_);
			auto const f2(score_(y2));

			boost::numeric::ublas::matrix<double> dfdv(shootfunc::NVAR, shootfunc::NVAR);
            
            // x1で用いる境界条件を変える
			auto const funcx1 = [&]
			{
                auto const sav = v1_;
                v1_ += delv1_;

                auto y = load1_(x1, v1_);
                integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y, x1, xf, dx_);
                auto const f = score_(y);

				// NVAR個の合致条件にある偏微分を数値的に計算
				for (auto i = 0U; i < shootfunc::NVAR; i++) {
					dfdv(i, 0) = (f[i] - f1[i]) / delv1_;
				}

				// 境界におけるパラメータを格納
				v1_ = sav;
			};

			cilk_spawn funcx1();
			
            // 次にx2で用いる境界条件を変える
			{	
				auto const sav = v2_;
				v2_ += delv2_;

				auto y = load2_(x2, v2_);
				integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y, x2, xf, - dx_);
				auto const f = score_(y);

				for (auto i = 0U; i < shootfunc::NVAR; i++)
					dfdv(i, 1) = (f2[i] - f[i]) / delv2_;
			
				v2_ = sav;
			};

			cilk_sync;

			shootfunc::dblasvector f(shootfunc::NVAR), ff(shootfunc::NVAR);
			for (auto i = 0U; i < shootfunc::NVAR; i++) {
				f[i] = f1[i] - f2[i];
				ff[i] = - f[i];
			}

			boost::numeric::ublas::permutation_matrix<> pm(dfdv.size1());
			boost::numeric::ublas::lu_factorize(dfdv, pm);              // 自由パラメータに対する増分を求める
			boost::numeric::ublas::lu_substitute(dfdv, pm, ff);

            v1_ += ff[0];                                               // x1の境界でのパラメータ値の増分

            v2_ += ff[1];                                               // x2の境界でのパラメータ値の増分

			y1 = load1_(x1, v1_);

			dvector res1;
			
            res1.reserve(boost::numeric_cast<std::size_t>((xf - x1) / dx_) + 2);
			
			auto const funcx1toxf = [&]
			{
				// 得られた条件でx1...dxまで微分方程式を解く
				integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y1, x1, dx_, dx_ - x1, [&res1](shootfunc::state_type const & y, double const x)
				{ res1.push_back(y[0]);	});									// x1...dxの結果を得る
				res1.pop_back();

				// 得られた条件でdx...xfまで微分方程式を解く
				integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y1, dx_, xf + x1, dx_, [&res1](shootfunc::state_type const & y, double const x)
				{ res1.push_back(y[0]); });									// dx...xf + x1の結果を得る
			};

			cilk_spawn funcx1toxf();

			// 得られた条件でx2...xfまで微分方程式を解く
			y2 = load2_(x2, v2_);								
			
            dvector res2;
			res2.reserve(boost::numeric_cast<std::size_t>((x2 - xf) / dx_) + 1);
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y2, x2, xf - x1, - dx_, [&res2](const shootfunc::state_type & y, const double x)
			{ res2.push_back(y[0]); });									// x2...xf - x1の結果を得る

			cilk_sync;

            return createResult(res1, res2, x1, xf);
		}

        shootf::result_type shootf::createResult(dvector const & res1, dvector const & res2, double x1, double xf) const
		{
			auto const size = boost::numeric_cast<std::vector<double>::size_type>(res1.size() + res2.size() - 1);

			dvector xp, xptmp, yp;
            xp.reserve(size);
            xptmp.reserve(size - 1);
            yp.reserve(size);

			auto const xfindex = boost::numeric_cast<std::size_t>(xf / dx_); 

			for (auto i = 0U; i < xfindex; i++) {
				xptmp.push_back(i ? static_cast<double>(i)* dx_ : x1);
			}

			for (auto i = xfindex + 1U; i < size; i++) {
				xptmp.push_back(static_cast<double>(i)* dx_);
			}

			yp.assign(res1.begin(), res1.end() - 1);
			auto const s = res2.size();
			for (auto i = 1U; i < s; i++) {
				yp.push_back(res2[s - i - 1]);
			}

			BOOST_ASSERT(xptmp.size() == yp.size());

			std::unique_ptr<gsl_interp_accel, decltype(utility::gsl_interp_accel_deleter)>
				const acc(gsl_interp_accel_alloc(), utility::gsl_interp_accel_deleter);
			std::unique_ptr<gsl_spline, decltype(utility::gsl_spline_deleter)>
				const spline(gsl_spline_alloc(gsl_interp_cspline, yp.size()), utility::gsl_spline_deleter);

			gsl_spline_init(spline.get(), xptmp.data(), yp.data(), xptmp.size());

			for (auto i = 0U; i < size; i++) {
				xp.push_back(i ? static_cast<double>(i) * dx_ : x1);
			}

			auto iter2 = yp.begin();
			std::advance(
                iter2,
                std::distance(
                    iter2,
                    boost::find_if(
                        yp,
				        [&, s](double x) { return std::fabs(x - res2[s - 2]) < EPS; })));

			yp.insert(iter2, gsl_spline_eval(spline.get(), xf, acc.get()));

			BOOST_ASSERT(xp.size() == yp.size());

			return std::make_pair(std::move(xp), std::move(yp));
		}
	}
}
