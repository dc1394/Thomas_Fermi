/*! \file makerhoenergy.cpp
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの実装

    Copyright © 2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "makerhoenergy.h"
#include <cmath>                                // for std::exp, std::pow
#include <iostream>                             // for std::cout
#include <utility>                              // for std::get
#include <boost/format.hpp>                     // for boost::format
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi

namespace thomasfermi {
	namespace makerhoen {
        // #region コンストラクタ

        MakeRhoEnergy::MakeRhoEnergy(std::int32_t n, MakeRhoEnergy::parameter_type const & pt, double Z) :
			alpha_(std::pow(128.0 / (9.0 * std::pow(boost::math::constants::pi<double>(), 2)) * Z, 1.0 / 3.0)),
			Z_(Z),
			b_(32.0 / (9.0 * std::pow(boost::math::constants::pi<double>(), 3)) * Z_ * Z_),
			xvec_(std::get<1>(pt)),
			dx_(xvec_[2] - xvec_[1]),
            fp_(nullptr, std::fclose),
			gl_(n),
            pbeta_(std::get<0>(pt)),
            size_(xvec_.size()),
			max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] / alpha_ / dx_)),
            y_prime_0_(std::get<2>(pt))
        {
            s_ = 4.0 * boost::math::constants::pi<double>() / (gl_.qgauss(
                [this](double x) { return std::sqrt(x) * y(x) * std::sqrt(y(x)); },
				xvec_.front(),
				xvec_.back()) * Z_);
        }

        // #endregion コンストラクタ

        // #region publicメンバ関数

		void MakeRhoEnergy::saveresult()
		{
			std::cout << boost::format("Energy = %.15f (Hartree)\n") % makeenergy();
			saverho("rho.csv");
			saverhoTilde("rhoTilde.csv");
			savey("y.csv");
		}

		// #endregion publicメンバ関数

		// #region privateメンバ関数
		
		double MakeRhoEnergy::exactrho(double r) const noexcept
		{
			return 4.0 * r * r * std::pow(Z_, 3) * std::exp(-2.0 * Z_ * r);
		}

		double MakeRhoEnergy::exactrhoTilde(double r) const noexcept
		{
			return 4.0 * std::pow(Z_, 3) * std::exp(-2.0 * Z_ * r);
		}

        double MakeRhoEnergy::makeenergy() const noexcept
        {
            return 3.0 / 7.0 * alpha_ * std::pow(Z_, 7.0 / 3.0) * y_prime_0_;
        }

        double MakeRhoEnergy::rho(double x) const noexcept
        {
			return s_ * b_ * std::pow(1.0 / alpha_, 2) * std::sqrt(x) * y(x) * std::sqrt(y(x));
        }

		double MakeRhoEnergy::rhoTilde(double x) const noexcept
		{
			return s_ * b_ * (y(x) / x) * std::sqrt(y(x) / x);
		}

        void MakeRhoEnergy::saverho(std::string const & filename)
        {
			fp_.reset(std::fopen(filename.c_str(), "w"));

            for (auto i = 1; i <= max_; i++) {
                auto const r = static_cast<double>(i) * dx_;
                std::fprintf(fp_.get(), "%.15f, %.15f, %.15f\n", r, rho(alpha_ * r), exactrho(r));
            }
        }

		void MakeRhoEnergy::saverhoTilde(std::string const & filename)
        {
			fp_.reset(std::fopen(filename.c_str(), "w"));

            for (auto i = 1; i <= max_; i++) {
				auto const r = static_cast<double>(i) * dx_;
				std::fprintf(fp_.get(), "%.15f, %.15f, %.15f\n", r, rhoTilde(alpha_ * r), exactrhoTilde(r));
            }
        }
                
        void MakeRhoEnergy::savey(std::string const & filename)
        {
			fp_.reset(std::fopen(filename.c_str(), "w"));

			for (auto const x : xvec_) {
				std::fprintf(fp_.get(), "%.15f, %.15f\n", x, y(x));
			}
        }

        double MakeRhoEnergy::y(double x) const
        {
			return std::pow(x * pbeta_->operator()<femall::Element::First>(x) * pbeta_->operator()<femall::Element::First>(x), 1.0 / 3.0);
        }

        // #endregion privateメンバ関数
	}
}
