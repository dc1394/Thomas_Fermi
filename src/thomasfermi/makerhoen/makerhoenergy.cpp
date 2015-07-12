/*! \file makerhoenergy.cpp
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "makerhoenergy.h"
#include <cmath>                                // for std::pow
#include <iostream>                             // for std::cout
#include <utility>                              // for std::get
#include <boost/cast.hpp>                       // for boost::numeric_cast
#include <boost/format.hpp>                     // for boost::format
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi

namespace thomasfermi {
	namespace makerhoen {
        // #region コンストラクタ

        MakeRhoEnergy::MakeRhoEnergy(std::int32_t n, MakeRhoEnergy::parameter_type const & pt, double Z) :
			alpha_(std::pow(128.0 / (9.0 * std::pow(boost::math::constants::pi<double>(), 2)) * Z, 1.0 / 3.0)),
			xvec_(std::get<1>(pt)),
			dx_(xvec_[2] - xvec_[1]),
            fp_(nullptr, fcloser),
			gl_(n),
            pbeta_(std::get<0>(pt)),
            size_(xvec_.size()),
			max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] / alpha_ / dx_)),
            v1_(std::get<2>(pt)),
            Z_(Z)
        {
            auto const func = myfunctional::make_functional(
                [this](double x) { return std::sqrt(x) * y(x) * std::sqrt(y(x)); });

            s_ = 1.0 / gl_.qgauss(
                func,
                xvec_.front(),
                xvec_.back());
        }

        // #endregion コンストラクタ

        // #region publicメンバ関数

		void MakeRhoEnergy::saveresult()
		{
			std::cout << boost::format("Energy = %.15f") % makeEnergy() << '\n';
			saverho("rho.csv");
			saverhoTilde("rhoTilde.csv");
			savey("y.csv");
		}

		// #endregion publicメンバ関数

		// #region privateメンバ関数

        double MakeRhoEnergy::makeEnergy() const
        {
            return 3.0 / 7.0 * alpha_ * Z_ * Z_ * v1_;
        }

        double MakeRhoEnergy::rho(double x) const
        {
            return alpha_ * s_ * Z_ * std::sqrt(x) * y(x) * std::sqrt(y(x));
        }

		double MakeRhoEnergy::rhoTilde(double x) const
		{
            auto const b = 32.0 / (9.0 * std::pow(boost::math::constants::pi<double>(), 3)) * Z_ * Z_;

			return b * s_ * (y(x) / x) * std::sqrt(y(x) / x);
		}

        void MakeRhoEnergy::saverho(std::string const & filename)
        {
			fp_.reset(std::fopen(filename.c_str(), "w"));

            for (auto i = 1; i <= max_; i++) {
                auto const r = static_cast<double>(i) * dx_;
                std::fprintf(fp_.get(), "%.15f, %.15f\n", r, rho(alpha_ * r));
            }
        }

		void MakeRhoEnergy::saverhoTilde(std::string const & filename)
        {
			fp_.reset(std::fopen(filename.c_str(), "w"));

            for (auto i = 1; i <= max_; i++) {
                auto const r = static_cast<double>(i) * dx_;
				std::fprintf(fp_.get(), "%.15f, %.15f\n", r, rhoTilde(alpha_ * r));
            }
        }
                
        void MakeRhoEnergy::savey(std::string const & filename)
        {
			fp_.reset(std::fopen(filename.c_str(), "w"));

			for (auto x : xvec_) {
				std::fprintf(fp_.get(), "%.15f, %.15f\n", x, y(x));
			}
        }

        double MakeRhoEnergy::y(double x) const
        {
			return std::pow(x * pbeta_->operator()<femall::Element::First>(x)* pbeta_->operator()<femall::Element::First>(x), 1.0 / 3.0);
        }

        // #endregion privateメンバ関数
	}
}
