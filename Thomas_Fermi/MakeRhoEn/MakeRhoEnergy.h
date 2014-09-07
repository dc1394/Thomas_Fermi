/*! \file MakeRhoEnergy.h
    \brief y(x)から電子密度とエネルギーを計算するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#ifndef _MAKERHOENERGY_H_
#define _MAKERHOENERGY_H_

#pragma once

#include "../Beta.h"
#include "../gausslegendre/Gauss_Legendre.h"
#include <cstdint>          // for std::int32_t
#include <fstream>          // for std::ofstream
#include <memory>           // for std::shared_ptr
#include <iomanip>          // for std::setprecision
#include <tuple>            // for std::tuple
#include <utility>          // for std::get
#include <boost/cast.hpp>   // for boost::numeric_cast

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
            // #region コンストラクタ・

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

            //! 原子のエネルギーを求める
            /*!
            \return 原子のエネルギーを求める
            */
            double makeEnergy() const;

            //! 関数ρ~(x)の値を返す
            /*!
            \param x xの値
            \return ρ~(x)の値
            */
            double rhoTilde(double x) const;

        public:
            //! 計算結果をファイルに出力する
            /*!
            \param n Gauss-Legendreの分点
            \param usecilk Intel® Cilk™ Plusを使うかどうか
            \param usesimd SIMDを使うかどうか
            */
            void saveResult();

        private:
            //! 関数ρ(x)の値をファイルに書き込む
            /*!
            \param filename 書き込むファイル名
            */
            void saverho(std::string const & filename);

            //! 関数ρ~(x)の値をファイルに書き込む
            /*!
            \param コピー元のオブジェクト（禁止）
            \return コピー元のオブジェクト
            */
            void saverhoTilde(std::string const & filename);

            //! 関数y(x)の値をファイルに書き込む
            /*!
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

            //! A static const variable.
            /*!
                エネルギーを計算するときの定数の値
            */
            static double const ALPHA;

            //! A private const variable.
            /*!
                x方向のメッシュ
            */
            std::vector<double> const xvec_;

            //! A private const variable.
            /*!
                x方向のメッシュの刻み幅
            */
            double const dx_;

            //! A private variable.
            /*!
                規格化のための定数
                s_ = 1.0 / (∫(0～∞)√x[y(x)]^(3/2)dx)
            */
            double s_;

            //! A private const variable.
            /*!
                ファイル出力するときのループの最大数
            */
            std::int32_t const max_;

            //! A private variable.
            /*!
                ファイル出力用のストリーム
            */
            std::ofstream ofs_;

            //! A private const variable.
            /*!
                ファイル出力用のストリーム
            */
            double const Z_;

            //! A private const variable.
            /*!
                Betaクラスのオブジェクトへのスマートポインタ
            */
            std::shared_ptr<const FEM_ALL::Beta> const pbeta_;

            //! A private const variable.
            /*!
                x方向のメッシュが格納されたstd::vectorのサイズ
            */
            std::size_t const size_;

            // #endregion メンバ変数

            // #region コンストラクタ

        public:
            //!  A constructor.
            /*!
            \param pt パラメータ
            \param Z 原子番号
            */
            MakeRhoEnergy(parameter_type const & pt, double Z);

            // #endregion コンストラクタ

            // #region operator=()

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
			:	xvec_(std::get<1>(pt)),
                dx_((xvec_[2] - xvec_[1]) * 2.0),
                max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] / MakeRhoEnergy::ALPHA / dx_)),
                pbeta_(std::get<0>(pt)),
				size_(xvec_.size()),
                Z_(Z)
		{
            gausslegendre::Gauss_Legendre gl(n);

            // 例外指定
            ofs_.exceptions(std::ios::badbit | std::ios::failbit);
            ofs_.setf(std::ios::fixed, std::ios::floatfield);
            ofs_ << std::setprecision(15);

            auto const func = myfunctional::make_functional(
                [&](double x) { return std::sqrt(x) * rhoTilde(x); });

            s_ = 1.0 / gl.qgauss(
                func,
                usesimd,
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
