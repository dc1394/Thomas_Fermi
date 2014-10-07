/*! \file MakeRhoEnergy.h
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/

#ifndef _MAKERHOENERGY_H_
#define _MAKERHOENERGY_H_

#pragma once

#include "../Beta.h"
#include "../gausslegendre/Gauss_Legendre.h"
#include <cstdint>                              // for std::int32_t
#include <fstream>                              // for std::ofstream
#include <memory>                               // for std::shared_ptr
#include <iomanip>                              // for std::setprecision
#include <tuple>                                // for std::tuple
#include <utility>                              // for std::get
#include <boost/cast.hpp>                       // for boost::numeric_cast
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi

namespace thomasfermi {
	namespace makerhoen {
        //! A class.
        /*!
            y(x)から電子密度とエネルギーを計算する
        */
		class MakeRhoEnergy final
		{
            // #region 型エイリアス

            typedef std::tuple<
                std::shared_ptr<const FEM_ALL::Beta> const,
                std::vector<double> const
                    > parameter_type;

            // #endregion 型エイリアス

        public:
            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
                \param n Gauss-Legendreの分点
                \param pt std::shared_ptr<Beta> constとstd::vector<double> constのタプル
                \param usesimd SIMDを使うかどうか
                \param Z 原子番号
            */
            MakeRhoEnergy(std::int32_t n, parameter_type const & pt, bool usesimd, double Z);

            //! A destructor.
            /*!
            */
            ~MakeRhoEnergy() {}

            // #endregion コンストラクタ・デストラクタ

        private:
            // #region メンバ関数

            //! A private member function (const).
            /*!
                原子のエネルギーを求める
                \return 原子のエネルギーを求める
            */
            double makeEnergy() const;

            //! A private member function (const).
            /*!
                xを引数にとり、関数ρ(x)の値を返す
                \param x xの値
                \return ρ(x)の値
            */
            double rho(double x) const;

            //! A private member function (const).
            /*!
                xを引数にとり、関数ρ~(x)の値を返す
                \param x xの値
                \return ρ~(x)の値
            */
            double rhoTilde(double x) const;

        public:
            //! A public member function.
            /*!
                計算結果をファイルに出力する
                \param n Gauss-Legendreの分点
                \param usecilk Intel® Cilk™ Plusを使うかどうか
                \param usesimd SIMDを使うかどうか
            */
            void saveResult();

        private:
            //! A private member function.
            /*!
                関数ρ(x)の値をファイルに書き込む
                \param filename 書き込むファイル名
            */
            void saverho(std::string const & filename);

            //! A private member function.
            /*!
                関数ρ~(x)の値をファイルに書き込む
                \param コピー元のオブジェクト（禁止）
                \return コピー元のオブジェクト
            */
            void saverhoTilde(std::string const & filename);

            //! A private member function.
            /*!
                関数y(x)の値をファイルに書き込む
                \param filename 書き込むファイル名
            */
            void savey(std::string const & filename);

            //! 関数y(x)の値を返す
            /*!
                \param x xの値
                \return y(x)の値
            */
            double y(double x) const;

            // #endregion メンバ関数

            // #region メンバ変数

            //! A private variable (constant).
            /*!
                エネルギーを計算するときの定数の値
            */
            double const alpha_;

            //! A private variable (constant).
            /*!
                x方向のメッシュが格納された動的配列
            */
            std::vector<double> const xvec_;

            //! A private variable (constant).
            /*!
                x方向のメッシュの刻み幅
            */
            double const dx_;

            //! A private variable (constant).
            /*!
                Gauss-Legendreの積分を行うオブジェクト
            */
            gausslegendre::Gauss_Legendre const gl_;

            //! A private variable.
            /*!
                規格化のための定数
                s_ = 1.0 / (∫(0～∞)√x[y(x)]^(3/2)dx)
            */
            double s_;

            //! A private variable (constant).
            /*!
                SIMDを使うかどうか
            */
            bool const usesimd_;

            //! A private variable (constant).
            /*!
                ファイル出力するときのループの最大数
            */
            std::int32_t const max_;

            //! A private variable.
            /*!
                ファイル出力用のストリーム
            */
            std::ofstream ofs_;

            //! A private variable (constant).
            /*!
                原子番号
            */
            double const Z_;

            //! A private variable (constant).
            /*!
                Betaクラスのオブジェクトへのスマートポインタ
            */
            std::shared_ptr<const FEM_ALL::Beta> const pbeta_;

            //! A private variable (constant).
            /*!
                x方向のメッシュが格納された動的配列のサイズ
            */
            std::size_t const size_;

            // #endregion メンバ変数

        private:
            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            MakeRhoEnergy() = delete;

            //! A private copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param コピー元のオブジェクト
            */
            MakeRhoEnergy(MakeRhoEnergy const &) = delete;

            //! operator=() (deleted).
            /*!
                operator=()（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            MakeRhoEnergy & operator=(MakeRhoEnergy const &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
		};

        inline MakeRhoEnergy::MakeRhoEnergy(std::int32_t n, parameter_type const & pt, bool usesimd, double Z)
            :   alpha_(std::pow(128.0 / (9.0 * power(boost::math::constants::pi<double>(), 2)) * Z, 1.0 / 3.0)),
                xvec_(std::get<1>(pt)),
                dx_((xvec_[2] - xvec_[1]) * 2.0),
                gl_(n),
                usesimd_(usesimd),
                max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] / alpha_ / dx_)),
                pbeta_(std::get<0>(pt)),
				size_(xvec_.size()),
                Z_(Z)
		{
            // 例外指定
            ofs_.exceptions(std::ios::badbit | std::ios::failbit);
            ofs_.setf(std::ios::fixed, std::ios::floatfield);
            ofs_ << std::setprecision(15);

            auto const func = myfunctional::make_functional(
                [&](double x) { return std::sqrt(x) * y(x) * std::sqrt(y(x)); });

            s_ = 1.0 / gl_.qgauss(
                func,
                usesimd_,
                xvec_[0],
                xvec_[size_ - 1]);
		}

        //! x^yを計算する（非メンバ関数）
        /*!
            \param x xの値    
            \param y yの値
            \return x^yの値
        */
		constexpr double power(double x, std::uint32_t y);
	}
}

#endif  // _MAKERHOENERGY_H_
