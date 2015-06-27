/*! \file foelementbeta.h
    \brief 一次要素でβ(x)を計算するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _FOELEMENTBETA_H_
#define _FOELEMENTBETA_H_

#pragma once

#include "beta.h"
#include <memory>   // for std::shared_ptr

namespace thomasfermi {
    namespace femall {
        //! A class.
        /*!
            一次要素でβ(x)を計算するクラス
        */
        class FoelementBeta final {
        public:
            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
                コンストラクタ
                \param pbeta Betaクラスのオブジェクトのスマートポインタ
            */
            FoelementBeta() = default;

            //! A destructor.
            /*!
                デフォルトデストラクタ
            */
            ~FoelementBeta() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ関数

            //! A public member function (const).
            /*!
                operator()(double x)の宣言
                一次要素でβ(x)を計算して返す
                \param x xの値
                \return 一次要素でのβ(x)の値
            */
            double operator()(Beta const & beta, double x) const;

            // #endregion メンバ関数

        private:
            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
            */
            FoelementBeta(FoelementBeta const &) = delete;

            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            FoelementBeta & operator=(FoelementBeta const &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };

        // #region メンバ関数

        inline double FoelementBeta::operator()(Beta const & beta, double x) const
        {
            auto klo = 0U;
            auto const max = static_cast<std::uint32_t>(beta.Size - 1);
            auto khi = max;

            // 表の中の正しい位置を二分探索で求める
            while (khi - klo > 1) {
                auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

                if ((*beta.Pxvec())[k] > x) {
                    khi = k;
                }
                else {
                    klo = k;
                }
            }

            // yvec[i] = f(xvec[i]), yvec[i + 1] = f(xvec[i + 1])の二点を通る直線を代入
            return ((*beta.Pyvec())[khi] - (*beta.Pyvec())[klo]) / ((*beta.Pxvec())[khi] - (*beta.Pxvec())[klo]) * (x - (*beta.Pxvec())[klo]) + (*beta.Pyvec())[klo];
        }
    }
}

#endif  // _FOELEMENTBETA_H_
