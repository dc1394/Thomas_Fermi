/*! \file shootf.cpp
    \brief 狙い撃ち法により、y(x)を求めるクラスの実装
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

#include "shootf.h"
#include <cmath>						// for std::fabs
#include <iterator>						// for std::advance, std::distance
#include <utility>                      // for std::make_pair
#include <boost/assert.hpp>				// for BOOST_ASSERT
#include <boost/range/algorithm.hpp>	// for boost::find_if			
#include <Eigen/Dense>
#include <Eigen/LU>						// for Eigen::FullPivLU
#include <gsl/gsl_spline.h>             // for gsl_interp_accel, gsl_interp_accel_free, gsl_spline, gsl_spline_free

#if _OPENMP >= 200805
    #include <omp.h>
#endif

namespace thomasfermi {
	namespace shoot {
        // #region コンストラクタ 

		shootf::shootf(double delv1, double delv2, double dx, double eps, loadfunctype const & load1, load2 const & l2,	scorefunctype const & score, double v1,	double v2)
	        :   delv1_(delv1),
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

        // #endregion コンストラクタ

        // #region publicメンバ関数

		void shootf::operator()(bool useomp, double x1, double x2, double xf, shootf::result_type & result)
		{
            using namespace boost::numeric::odeint;
            
			BOOST_ASSERT(x1 < xf);
			BOOST_ASSERT(x2 > xf);

			// 最良の仮の値v1_でx1からxfまで解いていく
			auto y1 = load1_(x1, v1_);
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y1, x1, xf, dx_);
			auto const f1(score_(y1));

			// 最良の仮の値v2_でx2からxfまで解いていく			
			auto y2 = load2_(x2, v2_);
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y2, x2, xf, - dx_);
			auto const f2(score_(y2));

			Eigen::MatrixXd dfdv(shootfunc::NVAR, shootfunc::NVAR);

            // x1で用いる境界条件を変える
			auto const funcx1 = [this](auto & dfdv, auto & f1, auto x1, auto xf)
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

            // 次にx2で用いる境界条件を変える    
            auto const functmp = [this](auto & dfdv, auto & f2, auto x2, auto xf)
            {
                using namespace boost::numeric::odeint;

				auto const sav = v2_;
				v2_ += delv2_;

				auto y = load2_(x2, v2_);
				integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y, x2, xf, - dx_);
				auto const f = score_(y);

				for (auto i = 0U; i < shootfunc::NVAR; i++) {
					dfdv(i, 1) = (f2[i] - f[i]) / delv2_;
				}
			
				v2_ = sav;
			};

			if (useomp) {
#if _OPENMP >= 200805
    #pragma omp task shared(dfdv)
#endif
                funcx1(dfdv, f1, x1, xf);
#if _OPENMP >= 200805
    #pragma omp task shared(dfdv)
#endif
                functmp(dfdv, f2, x2, xf);
#if _OPENMP >= 200805
    #pragma omp taskwait
#endif
			}
			else {
				funcx1(dfdv, f1, x1, xf);
                functmp(dfdv, f2, x2, xf);
			}
            
			Eigen::VectorXd f(shootfunc::NVAR), ff(shootfunc::NVAR);
			for (auto i = 0U; i < shootfunc::NVAR; i++) {
				f[i] = f1[i] - f2[i];
				ff[i] = - f[i];
			}

			Eigen::FullPivLU<Eigen::MatrixXd> lu(dfdv);
			ff = lu.solve(ff);
						
            v1_ += ff[0];                       // x1の境界でのパラメータ値の増分

            v2_ += ff[1];                       // x2の境界でのパラメータ値の増分

			y1 = load1_(x1, v1_);

			std::vector<double> res1, res2;
			if (useomp) {
#if _OPENMP >= 200805
    #pragma omp task shared(res1)
#endif
				res1 = solveodex1toxfpx1(x1, xf, y1);
#if _OPENMP >= 200805
    #pragma omp task shared(res2)
#endif
                res2 = solveodex2toxfmx1(x1, x2, xf);
#if _OPENMP >= 200805
    #pragma omp taskwait
#endif
			}
			else {
				res1 = solveodex1toxfpx1(x1, xf, y1);
                res2 = solveodex2toxfmx1(x1, x2, xf);
			}

            result = createResult(res1, res2, x1, xf);
		}

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        shootf::result_type shootf::createResult(std::vector<double> const & res1, std::vector<double> const & res2, double x1, double xf) const
		{
            auto const size = boost::numeric_cast<std::vector<double>::size_type>(res1.size() + res2.size() - 1);

			std::vector<double> xp, xptmp, yp;
            xp.reserve(size);
            xptmp.reserve(size - 1);
            yp.reserve(size);

			auto const xfindex = boost::numeric_cast<std::size_t>(xf / dx_); 

			for (auto i = 0U; i < xfindex; i++) {
				xptmp.push_back(i ? static_cast<double>(i) * dx_ : x1);
			}

			for (auto i = xfindex + 1U; i < size; i++) {
				xptmp.push_back(static_cast<double>(i) * dx_);
			}

			yp.assign(res1.begin(), res1.end() - 1);
			auto const s = res2.size();
			for (auto i = 1U; i < s; i++) {
				yp.push_back(res2[s - i - 1]);
			}

			BOOST_ASSERT(xptmp.size() == yp.size());

			std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>
				const acc(gsl_interp_accel_alloc(), gsl_interp_accel_free);
			std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)>
				const spline(gsl_spline_alloc(gsl_interp_cspline, yp.size()), gsl_spline_free);

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
				        [&res2, s](auto x) { return std::fabs(x - res2[s - 2]) < EPS; })));

			yp.insert(iter2, gsl_spline_eval(spline.get(), xf, acc.get()));

			BOOST_ASSERT(xp.size() == yp.size());

			return std::make_pair(std::move(xp), std::move(yp));
		}

        std::vector<double> shootf::solveodex2toxfmx1(double x1, double x2, double xf) const
        {
            using namespace boost::numeric::odeint;
            
			// 得られた条件でx2...xf - x1まで微分方程式を解く
			auto y2 = load2_(x2, v2_);								
			
            std::vector<double> res;
			res.reserve(boost::numeric_cast<std::size_t>((x2 - xf) / dx_) + 1);
                
            integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y2, x2, xf - x1, - dx_, [&res](auto const & y, auto const)
			{ res.push_back(y[0]); });         // x2...xf - x1の結果を得る
			
            return res;
        }

        std::vector<double> shootf::solveodex1toxfpx1(double x1, double xf, shootfunc::state_type & y1) const
        {
            using namespace boost::numeric::odeint;

            std::vector<double> res;
            res.reserve(boost::numeric_cast<std::size_t>((xf - x1) / dx_) + 2);
            
            // 得られた条件でx1...dxまで微分方程式を解く
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y1, x1, dx_, dx_ - x1, [&res](auto const & y, auto const)
			{ res.push_back(y[0]); });     // x1...dxの結果を得る
			res.pop_back();
            
			// 得られた条件でdx...xf + x1まで微分方程式を解く
			integrate_const(stepper_type(eps_, eps_), shootfunc::rhs, y1, dx_, xf + x1, dx_, [&res](auto const & y, auto const)
			{ res.push_back(y[0]); });     // dx...xf + x1の結果を得る
            
            //std::printf("%d\n", static_cast<int>(res.size()));
            return res;
        }

        // #endregion privateメンバ関数
	}
}

