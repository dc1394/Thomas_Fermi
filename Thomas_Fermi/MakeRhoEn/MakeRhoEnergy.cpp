#include "MakeRhoEnergy.h"
#include "../Gauss_Legendre.h"
#include <iostream>
#include <cmath>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

namespace Thomas_Fermi {
	namespace MakeRhoEn {
		const double MakeRhoEnergy::ALPHA = std::pow(128.0 / (9.0 * power(boost::math::constants::pi<double>(), 2)), 1.0 / 3.0);

		double MakeRhoEnergy::y(double x) const
		{
			return std::pow(x * (*pbeta_)(x) * (*pbeta_)(x), 1.0 / 3.0);
		}

		double MakeRhoEnergy::rhofunc(double x) const
		{
			static const double fixednum = 32.0 / (9.0 * power(boost::math::constants::pi<double>(), 3));
			
			return fixednum * (y(x) / x) * std::sqrt(y(x) / x);
		}

		double MakeRhoEnergy::makeEnergy() const
		{
			//return 3.0 / 7.0 * ALPHA * v1_;
		}

        void MakeRhoEnergy::savey(const std::string & filename)
		{
			ofs_.open(filename);

			for (std::uint32_t i = 0; i < size_; i++)
				ofs_ << xvec_[i] << ' ' << y(xvec_[i]) << '\n';

			ofs_.close();
		}

		void MakeRhoEnergy::saverho(const std::string & filename)
		{
			ofs_.open(filename);

			for (std::int32_t i = 1; i <= max_; i++) {
				const double r = static_cast<double>(i) * dx_;
				ofs_ << r << boost::format(",%.15f,%.15f\n") % (rhofunc(MakeRhoEnergy::ALPHA * r) / a);
			}

			ofs_.close();
		}

		void MakeRhoEnergy::saveresult3()
		{
			ofs_.open("result3.csv");

			for (std::int32_t i = 1; i <= max_; i++) {
				const double r = static_cast<double>(i) * dx_;
				ofs << r << boost::format(",%.15f,%.15f\n") % (r * r * rhofunc(MakeRhoEnergy::ALPHA * r) / a)
															% (4.0 * r * r * std::exp(- 2.0 * r));
			}

			ofs_.close();
		}

        void MakeRhoEnergy::saveresult(std::size_t n, bool usecilk, bool usesimd, double Z)
		{
            FEM_ALL::Gauss_Legendre gl(n, usecilk);

            const double a = 128.0 * Z * Z / (9.0 * power(boost::math::constants::pi<double>(), 2));

            k_ = a / gl.qgauss(
                xvec_[0],
                xvec_[size_ - 1],
                usesimd,
                [&](double x) { return x * x * rhofunc(x); });

			std::cout << boost::format("Energy = %.15f * Z^(7/3) (Hartree)\n") % makeEnergy();

			ofs_.open("Energy.txt");
			
			ofs << boost::format("y'(0) = %.15f\n") % v1_;
			ofs << boost::format("Energy = %.15f * Z^(7/3) (Hartree)\n") % makeEnergy();
			
			ofs.close();

			saveresult1();
			saveresult2(n, useAVXorSSE, usecilk);
			saveresult3();
		}

		constexpr double power(double x, std::uint32_t y)
		{
			return y == 0 ? 1 : power(x, y - 1) * x;
		}
	}
}
