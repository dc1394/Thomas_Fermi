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
#ifdef __GXX_EXPERIMENTAL_CXX0X__
			static constexpr double fixednum = 32.0 / (9.0 * power(boost::math::constants::pi<double>(), 3));
#else
			static const double fixednum = 32.0 / (9.0 * power(boost::math::constants::pi<double>(), 3));
#endif			
			return fixednum * (y(x) / x) * std::sqrt(y(x) / x);
		}

		double MakeRhoEnergy::makeEnergy() const
		{
			return 3.0 / 7.0 * ALPHA * v1_;
		}

		void MakeRhoEnergy::saveresult1()
		{
			ofs.open("result1.csv");

			for (std::uint32_t i = 0; i < size_; i++)
				ofs << xvec_[i] << boost::format(",%.15f\n") % y(xvec_[i]);

			ofs.close();
		}

		void MakeRhoEnergy::saveresult2(std::size_t n, bool useAVXorSSE, bool usecilk)
		{
			ofs.open("result2.csv");

			FEM_ALL::Gauss_Legendre gl(n, usecilk);
			a = gl.qgauss(xvec_[0], xvec_[size_ - 1], useAVXorSSE,
					      [&](double x) { return x * x * rhofunc(x); }) / 
							power(MakeRhoEnergy::ALPHA, 3);

			for (std::int32_t i = 1; i <= max_; i++) {
				const double r = static_cast<double>(i) * dx_;
				ofs << r << boost::format(",%.15f,%.15f\n") % (rhofunc(MakeRhoEnergy::ALPHA * r) / a)
															% (4.0 * std::exp(- 2.0 * r));
			}

			ofs.close();
		}

		void MakeRhoEnergy::saveresult3()
		{
			ofs.open("result3.csv");

			for (std::int32_t i = 1; i <= max_; i++) {
				const double r = static_cast<double>(i) * dx_;
				ofs << r << boost::format(",%.15f,%.15f\n") % (r * r * rhofunc(MakeRhoEnergy::ALPHA * r) / a)
															% (4.0 * r * r * std::exp(- 2.0 * r));
			}

			ofs.close();
		}

		void MakeRhoEnergy::saveresult(std::size_t n, bool useAVXorSSE, bool usecilk)
		{
			// 例外指定
			ofs.exceptions(std::ios::badbit | std::ios::failbit);

			std::cout << boost::format("y'(0) = %.15f\n") % v1_;
			std::cout << boost::format("Energy = %.15f * Z^(7/3) (Hartree)\n") % makeEnergy();

			ofs.open("Energy.txt");
			
			ofs << boost::format("y'(0) = %.15f\n") % v1_;
			ofs << boost::format("Energy = %.15f * Z^(7/3) (Hartree)\n") % makeEnergy();
			
			ofs.close();

			saveresult1();
			saveresult2(n, useAVXorSSE, usecilk);
			saveresult3();
		}

#ifdef __GXX_EXPERIMENTAL_CXX0X__
		constexpr
#endif
		double power(double x, std::uint32_t y)
		{
			return y == 0 ? 1 : power(x, y - 1) * x;
		}
	}
}
