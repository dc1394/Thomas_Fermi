/*! \file MakeRhoEnergy.cpp
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "MakeRhoEnergy.h"
#include "../Gauss_Legendre.h"
#include <cmath>                                // for std::pow
#include <iostream>                             // for std::cout
#include <utility>                              // for std::get
#include <boost/assert.hpp>                     // for boost::assert
#include <boost/cast.hpp>                       // for boost::numeric_cast
#include <boost/format.hpp>                     // for boost::format
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi

namespace thomasfermi {
	namespace makerhoen {
        // #region コンストラクタ

        MakeRhoEnergy::MakeRhoEnergy(std::int32_t n, MakeRhoEnergy::parameter_type const & pt, bool usesimd, double Z)
            : alpha_(std::pow(128.0 / (9.0 * std::pow(boost::math::constants::pi<double>(), 2)) * Z, 1.0 / 3.0)),
            xvec_(std::get<1>(pt)),
            dx_((xvec_[2] - xvec_[1]) * 2.0),
            gl_(n),
            usesimd_(usesimd),
            max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] / alpha_ / dx_)),
            pbeta_(std::get<0>(pt)),
            size_(xvec_.size()),
            Z_(Z)
        {
            auto const func = myfunctional::make_functional(
                [&](double x) { return std::sqrt(x) * y(x) * std::sqrt(y(x)); });

            s_ = 1.0 / gl_.qgauss(
                func,
                usesimd_,
                xvec_[0],
                xvec_[size_ - 1]);
        }

        // #endregion コンストラクタ

        // #region privateメンバ関数

        // #region publicメンバ関数


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

        // #endregion privateメンバ関数

        


	}
}
