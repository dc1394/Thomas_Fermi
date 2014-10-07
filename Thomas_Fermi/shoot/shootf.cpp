#include "shootf.h"
#include <cmath>
#include <cstdint>
#include <iterator>
#include <boost/assert.hpp>
#include <boost/cast.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/range/algorithm.hpp>

namespace thomasfermi {
	namespace shoot {
		shootf::result_type shootf::operator()(double x1, double x2, double xf)
		{
			BOOST_ASSERT(x1 < xf);
			BOOST_ASSERT(x2 > xf);

			using namespace boost::numeric::odeint;

			bulirsch_stoer<shootfunc::state_type> stepper(eps_, eps_);

			auto y1 = load1_(v1_, x1);					                // 最良の仮の値v1_でx1からxfまで解いていく
			integrate_const(stepper, shootfunc::rhs, y1, x1, xf, dx_);
			auto const f1(score_(y1));
		
			auto y2 = load2_(v2_, x2);					                // 最良の仮の値v2_でx2からxfまで解いていく			
			integrate_const(stepper, shootfunc::rhs, y2, x2, xf, - dx_);
			auto const f2(score_(y2));

			boost::numeric::ublas::matrix<double> dfdv(shootfunc::NVAR, shootfunc::NVAR);
            
            // x1で用いる境界条件を変える
            {
                auto const sav = v1_;
                v1_ += delv1_;

                auto const y = load1_(v1_, x1);
                integrate_const(stepper, shootfunc::rhs, y, x1, xf, dx_);
                auto const f = score_(y);

                for (std::size_t i = 0; i < shootfunc::NVAR; i++)		// NVAR個の合致条件にある偏微分を数値的に計算
                    dfdv(i, 0) = (f[i] - f1[i]) / delv1_;

                v1_ = sav;											    // 境界におけるパラメータを格納
            }

            // 次にx2で用いる境界条件を変える
			{	
				auto const sav = v2_;
				v2_ += delv2_;

				auto const y = load2_(v2_, x2);
				integrate_const(stepper, shootfunc::rhs, y, x2, xf, - dx_);
				auto const f = score_(y);

				for (std::size_t i = 0; i < shootfunc::NVAR; i++)
					dfdv(i, 1) = (f2[i] - f[i]) / delv2_;
			
				v2_ = sav;
			}

			shootfunc::dblasvector f(shootfunc::NVAR), ff(shootfunc::NVAR);
			for (std::size_t i = 0; i < shootfunc::NVAR; i++) {
				f[i] = f1[i] - f2[i];
				ff[i] = - f[i];
			}

			boost::numeric::ublas::permutation_matrix<> pm(dfdv.size1());
			boost::numeric::ublas::lu_factorize(dfdv, pm);              // 自由パラメータに対する増分を求める
			boost::numeric::ublas::lu_substitute(dfdv, pm, ff);

            v1_ += ff[0];                                               // x1の境界でのパラメータ値の増分

            v2_ += ff[1];                                               // x2の境界でのパラメータ値の増分

			y1 = load1_(v1_, x1);

			dvector res1;
			
            res1.reserve(boost::numeric_cast<std::size_t>((xf - x1) / dx_) + 2);
			
            // 得られた条件でx1...dxまで微分方程式を解く
			integrate_const(stepper, shootfunc::rhs, y1, x1, dx_, dx_ - x1, [&res1](shootfunc::state_type const & y, double const x)
			{ res1.push_back(y[0]); });									// x1...dxの結果を得る
			res1.pop_back();

			// 得られた条件でdx...xfまで微分方程式を解く
			integrate_const(stepper, shootfunc::rhs, y1, dx_, xf + x1, dx_, [&res1](shootfunc::state_type const & y, double const x)
			{ res1.push_back(y[0]); });									// dx...xf + x1の結果を得る
			
			// 得られた条件でx2...xfまで微分方程式を解く
			y2 = load2_(v2_, x2);								
			
            dvector res2;
			res2.reserve(boost::numeric_cast<std::size_t>((x2 - xf) / dx_) + 1);
			integrate_const(stepper, shootfunc::rhs, y2, x2, xf - x1, - dx_, [&res2](const shootfunc::state_type & y, const double x)
			{ res2.push_back(y[0]); });									// x2...xf - x1の結果を得る

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

			for (std::size_t i = 0; i < xfindex; i++)
				xptmp.push_back(i ? static_cast<double>(i) * dx_ : x1);

			for (std::size_t i = xfindex + 1; i < size; i++)
				xptmp.push_back(static_cast<double>(i) * dx_);

			yp.assign(res1.begin(), res1.end() - 1);

			auto const s = res2.size();
			
            for (std::size_t i = 1; i < s; i++)
				yp.push_back(res2[s - i - 1]);

			BOOST_ASSERT(xptmp.size() == yp.size());

            alglib::spline1dinterpolant spline;

            alglib::real_1d_array x, y;

            x.setcontent(xptmp.size(), xptmp.data());
            y.setcontent(yp.size(), yp.data());

            alglib::spline1dbuildakima(x, y, spline);
			
			for (std::size_t i = 0; i < size; i++)
				xp.push_back(i ? static_cast<double>(i) * dx_ : x1);

			auto const iter2 = yp.begin();
			std::advance(
                iter2,
                std::distance(
                    iter2,
                    boost::find_if(
                        yp,
				        [&, s](double x) { return std::fabs(x - res2[s - 2]) < EPS; })));

            yp.insert(iter2, alglib::spline1dcalc(spline, xf));

			BOOST_ASSERT(xp.size() == yp.size());

			return std::make_tuple(std::move(xp), std::move(yp));
		}
	}
}
