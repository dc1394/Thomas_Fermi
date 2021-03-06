/*! \file MakeRhoEnergy.h
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _MAKERHOENERGY_H_
#define _MAKERHOENERGY_H_

#pragma once

#include "../beta.h"
#include "../gausslegendre/gausslegendre.h"
#include <cstdint>                          // for std::int32_t
#include <cstdio>                           // for FILE, std::fclose
#include <memory>                           // for std::shared_ptr
#include <tuple>                            // for std::tuple

namespace thomasfermi {
    namespace makerhoen {
        //! A class.
        /*!
            y(x)から電子密度とエネルギーを計算する
        */
        class MakeRhoEnergy final {
            // #region 型エイリアス

            using parameter_type = std::tuple<std::shared_ptr<femall::Beta>, std::vector<double>, double const>;

            // #endregion 型エイリアス

        public:
            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
                唯一のコンストラクタ
                \param n Gauss-Legendreの分点
                \param pt std::vector<double>、std::shared_ptr<Beta>、doubleのstd::tuple
                \param Z 原子番号
            */
            MakeRhoEnergy(std::int32_t n, parameter_type const & pt, double Z);

            //! A default destructor.
            /*!
                デフォルトデストラクタ
            */
            ~MakeRhoEnergy() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region publicメンバ関数

            //! A public member function.
            /*!
                計算結果をファイルに出力する
            */
            void saveresult();

            // #endregion publicメンバ関数

        private:
            // #region privateメンバ関数

            //! A private member function (const).
            /*!
                厳密な電子密度の関数
                \param r 原点からの距離（原子単位）
                \return 厳密な電子密度
            */
            double exactrho(double r) const noexcept;

            //! A private member function (const).
            /*!
                厳密な電子密度の関数（4πr ** 2で割っている）
                \param r 原点からの距離（原子単位）
                \return 厳密な電子密度の関数（4πr ** 2で割っている）
            */
            double exactrhoTilde(double r) const noexcept;

            //! A private member function (const).
            /*!
                原子のエネルギーを求める
                \return 原子のエネルギーを求める
            */
            double makeenergy() const noexcept;

            //! A private member function (const).
            /*!
                xを引数にとり、関数ρ(x)の値を返す
                \param x xの値
                \return ρ(x)の値
            */
            double rho(double x) const noexcept;

            //! A private member function (const).
            /*!
                xを引数にとり、関数ρ~(x)の値を返す
                \param x xの値
                \return ρ~(x)の値
            */
            double rhoTilde(double x) const noexcept;

            //! A private member function.
            /*!
                関数ρ(x)の値をファイルに書き込む
                \param filename 書き込むファイル名
            */
            void saverho(std::string const & filename);

            //! A private member function.
            /*!
                関数ρ~(x)の値をファイルに書き込む
                \param filename コピー元のオブジェクト（禁止）
                \return コピー元のオブジェクト
            */
            void saverhoTilde(std::string const & filename);

            //! A private member function.
            /*!
                関数y(x)の値をファイルに書き込む
                \param filename 書き込むファイル名
            */
            void savey(std::string const & filename);

            //! A private member function.
            /*!
                関数y(x)の値を返す
                \param x xの値
                \return y(x)の値
            */
            double y(double x) const;

            // #endregion privateメンバ関数

            // #region メンバ変数

            //! A private variable (constant).
            /*!
                α = [128 / (9π ** 2)]^(1 / 3) * Z^(1 / 3)
            */
            double const alpha_;

            //! A private variable (constant).
            /*!
                原子番号
            */
            double const Z_;

            //! A private variable (constant).
            /*!
                b = 32 / (9π ** 3) * Z ** 2
            */
            double const b_;

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
                ファイルポインタ
            */
            std::unique_ptr<FILE, decltype(&std::fclose)> fp_;

            //! A private variable (constant).
            /*!
                Gauss-Legendreの積分を行うオブジェクト
            */
            gausslegendre::Gauss_Legendre const gl_;

            //! A private variable (constant).
            /*!
                Betaクラスのオブジェクトへのスマートポインタ
            */
            std::shared_ptr<femall::Beta> const pbeta_;

            //! A private variable.
            /*!
                規格化のための定数
                s_ = 1.0 / (∫(0～∞)√x[y(x)]^(3/2)dx)
            */
            double s_;

            //! A private variable (constant).
            /*!
                x方向のメッシュが格納された動的配列のサイズ
            */
            std::size_t const size_;

            //! A private variable (constant).
            /*!
                ファイル出力するときのループの最大数
            */
            std::int32_t const max_;

            //! A private member variable.
            /*!
                原点に近いxにおけるyの微分値
            */
            double const y_prime_0_;

            // #endregion メンバ変数

        public:
            // #region 禁止されたコンストラクタ・メンバ関数

            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            MakeRhoEnergy() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト
            */
            MakeRhoEnergy(MakeRhoEnergy const & dummy) = delete;

            //! operator=() (deleted).
            /*!
                operator=()（禁止）
                \param dummy コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            MakeRhoEnergy & operator=(MakeRhoEnergy const & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }
}

#endif  // _MAKERHOENERGY_H_
