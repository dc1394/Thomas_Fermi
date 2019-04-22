/*! \file beta.h
    \brief β(x)を計算するクラスの宣言
    Copyright © 2014-2019 @dc1394 All Rights Reserved.
   
    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BETA_H_
#define _BETA_H_

#pragma once

#include "element.h"
#include "utility/property.h"
#include <cstdint>              // for std::uint32_t
#include <memory>               // for std::shared_ptr
#include <vector>               // for std::vector

namespace thomasfermi {
    namespace femall {
        //! A class.
        /*!
            β(x)を計算するクラス
        */
        class Beta final {
            // #region コンストラクタ・デストラクタ

        public:
            //! A constructor.
            /*!
                コンストラクタ
                \param xvec x方向のメッシュ
                \param yvec y方向のメッシュ
            */
            Beta(std::vector<double> const & xvec, std::vector<double> const & yvec);
            
            //! A destructor.
            /*!
                デフォルトデストラクタ
            */
            ~Beta() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ関数

            template <Element E>
            //! A public member function (const).
            /*!
                一次要素でβ(x)を計算して返す
                \param x xの値
                \return β(x)の値
            */
            double operator()(double x) const;

            template <Element E>
            //! A public member function (const).
            /*!
                一次要素でβ(x)/dxを計算して返す
                \param x xの値
                \return dβ(x)/dxの値
            */
            double dxdbeta(double x) const;

            // #endregion メンバ関数

        private:
            // #region メンバ変数

            //!  A private member variable (constant).
            /*!
                メッシュが格納された動的配列のサイズ
            */
            std::size_t const size_;
            
            //!  A private member variable (constant).
            /*!
                x方向のメッシュが格納された動的配列
            */
            std::vector<double> const xvec_;

            //!  A private member variable (constant).
            /*!
                y方向のメッシュが格納された動的配列
            */
            std::vector<double> const yvec_;

            // #endregion メンバ変数

            // #region 禁止されたコンストラクタ・メンバ関数

        public:
            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            Beta() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト（未使用）
            */
            Beta(Beta const & dummy) = delete;

            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            Beta & operator=(Beta const &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };

        // #region コンストラクタ

        inline Beta::Beta(std::vector<double> const & xvec, std::vector<double> const & yvec) :
            size_(xvec.size()),
            xvec_(xvec),
            yvec_(yvec)
        {
        }

        // #endregion コンストラクタ

        // #region メンバ関数

        template <>
        inline double Beta::operator()<Element::First>(double x) const
        {
            auto klo = 0U;
            auto const max = static_cast<std::uint32_t>(size_ - 1);
            auto khi = max;

            // 表の中の正しい位置を二分探索で求める
            while (khi - klo > 1) {
                auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

                if (xvec_[k] > x) {
                    khi = k;
                }
                else {
                    klo = k;
                }
            }

            // yvec_[i] = f(xvec_[i]), yvec_[i + 1] = f(xvec_[i + 1])の二点を通る直線を代入
            return (yvec_[khi] - yvec_[klo]) / (xvec_[khi] - xvec_[klo]) * (x - xvec_[klo]) + yvec_[klo];
        }
               
        template <>
        inline double Beta::operator()<Element::Second>(double x) const
        {
            auto klo = 0U;
            auto const max = static_cast<std::uint32_t>(size_ - 1);
            auto khi = max;

            // 表の中の正しい位置を二分探索で求める
            while (khi - klo > 1) {
                auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

                if (xvec_[k] > x) {
                    khi = k;
                }
                else {
                    klo = k;
                }
            }
            
            // yvec_[i] = f(xvec_[i]), yvec_[i + 1] = f(xvec_[i + 1]), yvec_[i + 2] = f(xvec_[i + 2])の三点を通る放物線を生成
            
            // もし、配列の外にはみ出るときは
            if (khi >= max - 1) {
                // 一つ前の値を使う
                khi--;
                klo--;
            }

            auto const x2mx1 = xvec_[khi + 1] - xvec_[khi];
            auto const x0mx2 = xvec_[klo] - xvec_[khi + 1];
            auto const x1mx0 = xvec_[khi] - xvec_[klo];

            auto const denom = x2mx1 * x0mx2 * x1mx0;

            auto const a = -(x2mx1 * yvec_[klo] + x0mx2 * yvec_[khi] + x1mx0 * yvec_[khi + 1]);
            
            auto b = x2mx1 * (xvec_[khi + 1] + xvec_[khi]) * yvec_[klo];
            b += x0mx2 * (xvec_[klo] + xvec_[khi + 1]) * yvec_[khi];
            b += x1mx0 * (xvec_[khi] + xvec_[klo]) * yvec_[khi + 1];

            auto c = x2mx1 * xvec_[khi] * xvec_[khi + 1] * yvec_[klo];
            c += x0mx2 * xvec_[klo] * xvec_[khi + 1] * yvec_[khi];
            c += x1mx0 * xvec_[klo] * xvec_[khi] * yvec_[khi + 1];
            c *= -1.0;

            return ((a * x + b) * x + c) / denom;
        }

        template <>
        inline double Beta::dxdbeta<Element::First>(double x) const
        {
            auto klo = 0U;
            auto const max = static_cast<std::uint32_t>(size_ - 1);
            auto khi = max;

            // 表の中の正しい位置を二分探索で求める
            while (khi - klo > 1) {
                auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

                if (xvec_[k] > x) {
                    khi = k;
                }
                else {
                    klo = k;
                }
            }

            // yvec_[i] = f(xvec_[i]), yvec_[i + 1] = f(xvec_[i + 1])の二点を通る直線の傾き
            return (yvec_[khi] - yvec_[klo]) / (xvec_[khi] - xvec_[klo]);
        }

        // #endregion メンバ関数
    }
}

#endif  // _BETA_H_

