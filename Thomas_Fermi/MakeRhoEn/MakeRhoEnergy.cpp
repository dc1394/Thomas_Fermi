#include "MakeRhoEnergy.h"
#include "../Gauss_Legendre.h"

#include <cmath>                                // for std::pow
#include <iostream>                             // for std::cout
#include <boost/assert.hpp>                     // for boost::assert
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi

namespace Thomas_Fermi {
	namespace MakeRhoEn {
		const double MakeRhoEnergy::ALPHA = std::pow(128.0 / (9.0 * power(boost::math::constants::pi<double>(), 2)), 1.0 / 3.0);

        double MakeRhoEnergy::makeEnergy() const
        {
            //return 3.0 / 7.0 * ALPHA * v1_;
        }

		double MakeRhoEnergy::rho(double x) const
		{
            const double a = 128.0  * Z_ * Z_ / (9.0 * power(boost::math::constants::pi<double>(), 2));
			
			return a * s_ * (y(x) / x) * std::sqrt(y(x) / x);
		}

        void MakeRhoEnergy::saveresult(std::size_t n, bool usecilk, bool usesimd)
		{
            FEM_ALL::Gauss_Legendre gl(n, usecilk);

            s_ = 1.0 / gl.qgauss(
                xvec_[0],
                xvec_[size_ - 1],
                usesimd,
                [&](double x) { return std::sqrt(x) * rho(x); });

			std::cout << boost::format("Energy = %.15f * Z^(7/3) (Hartree)\n") % makeEnergy();

			ofs_.open("Energy.txt");
			
			ofs_ << boost::format("Energy = %.15f * Z^(7/3) (Hartree)\n") % makeEnergy();
			
			ofs_.close();

			saverho("rho.txt");
			savey("y_x.txt");
			saverhoTilde("rhoTilde.txt");
		}

        void MakeRhoEnergy::saverho(const std::string & filename)
        {
            ofs_.open(filename);

            for (std::int32_t i = 1; i <= max_; i++) {
                const double r = static_cast<double>(i)* dx_;
                ofs_ << r << ' ' << rho(MakeRhoEnergy::ALPHA * r);
            }

            ofs_.close();
        }

        void MakeRhoEnergy::saverhoTilde(const std::string & filename)
        {
            ofs_.open(filename);

            for (std::int32_t i = 1; i <= max_; i++) {
                const double r = static_cast<double>(i)* dx_;
                ofs_ << r << ' ' << r * r * rho(MakeRhoEnergy::ALPHA * r) << '\n';
            }

            ofs_.close();
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
            return std::pow(x * (*pbeta_)(x)* (*pbeta_)(x), 1.0 / 3.0);
        }

		constexpr double power(double x, std::uint32_t y)
		{
			return y == 0 ? 1 : power(x, y - 1) * x;
		}
	}
}
