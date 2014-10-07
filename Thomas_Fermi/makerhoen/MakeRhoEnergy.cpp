/*! \file MakeRhoEnergy.cpp
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/

#include "MakeRhoEnergy.h"
#include "../Gauss_Legendre.h"

#include <cmath>                                // for std::pow
#include <iostream>                             // for std::cout
#include <boost/assert.hpp>                     // for boost::assert
#include <boost/format.hpp>                     // for boost::format

namespace thomasfermi {
	namespace makerhoen {
        double MakeRhoEnergy::makeEnergy() const
        {
            auto const func = myfunctional::make_functional(
                [&](double x) { return (*pbeta_)(x); });

            auto const sum = gl_.qgauss(
                func,
                usesimd_,
                xvec_[0],
                xvec_[size_ - 1]);

            return -3.0 / 7.0 * alpha_ * Z_ * Z_ * sum;
        }

        double MakeRhoEnergy::rho(double x) const
        {
            return alpha_ * Z_ * s_ * std::sqrt(x) * y(x) * std::sqrt(y(x));
        }

		double MakeRhoEnergy::rhoTilde(double x) const
		{
            auto const b = 32.0 / (9.0 * power(boost::math::constants::pi<double>(), 3)) * Z_ * Z_;

			return b * s_ * (y(x) / x) * std::sqrt(y(x) / x);
		}

        void MakeRhoEnergy::saverho(const std::string & filename)
        {
            ofs_.open(filename);

            for (std::int32_t i = 1; i <= max_; i++) {
                auto const r = static_cast<double>(i) * dx_;
                ofs_ << r << ' ' << rho(alpha_ * r);
            }

            ofs_.close();
        }

        void MakeRhoEnergy::saverhoTilde(const std::string & filename)
        {
            ofs_.open(filename);

            for (std::int32_t i = 1; i <= max_; i++) {
                auto const r = static_cast<double>(i) * dx_;
                ofs_ << r << ' ' << rhoTilde(alpha_ * r) << '\n';
            }

            ofs_.close();
        }

        void MakeRhoEnergy::saveResult()
		{
			std::cout << boost::format("Energy = %.15f (Hartree)\n") % makeEnergy();

			ofs_.open("Energy.txt");
			
			ofs_ << boost::format("Energy = %.15f (Hartree)\n") % makeEnergy();
			
			ofs_.close();

			saverho("rho.txt");
            saverhoTilde("rhoTilde.txt");
			savey("y_x.txt");
		}
        
        void MakeRhoEnergy::savey(const std::string & filename)
        {
            ofs_.open(filename);

            for (auto x : xvec_)
                ofs_ << x << ' ' << y(x) << '\n';

            ofs_.close();
        }

        double MakeRhoEnergy::y(double x) const
        {
            return std::pow(x * (*pbeta_)(x) * (*pbeta_)(x), 1.0 / 3.0);
        }

		constexpr double power(double x, std::uint32_t y)
		{
			return y == 0 ? 1 : power(x, y - 1) * x;
		}
	}
}
