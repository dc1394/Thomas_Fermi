/*! \file shootfunc.h
    \brief y(x)の初期関数y0(x)の、原点に近い点xにおける関数値および微分値と、
           適合点xfにおける関数値および微分値を求めるクラスの宣言
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

#ifndef _SHOOTFUNC_H_
#define _SHOOTFUNC_H_

#pragma once

#include <array>        // for std::array
#include <Eigen/Core>   // Eigen::VectorXd

namespace thomasfermi {
    namespace shoot {
        //! A class.
        /*!
            y(x)の初期関数y0(x)の、原点に近い点xにおける関数値および微分値と、
            適合点xfにおける関数値および微分値を求めるクラス
        */
        class shootfunc final {
        public:
            // #region メンバ定数

            //! A public member variable (constant expression).
            /*!
                二元連立常微分方程式を表す定数
            */
            static std::size_t constexpr NVAR = 2;

            // #endregion メンバ定数

            // #region 型エイリアス

            using state_type = std::array<double, NVAR>;

            // #region 型エイリアス

            // #region メンバ関数

            //! A public static member function.
            /*!
                y0(x)の原点に近い点xにおける関数値および微分値を求める
                \param x1 xの値
                \param v1 y'0(0)
                \return y0(x)の原点に近い点xにおける関数値および微分値
            */
            static shootfunc::state_type load1(double x1, double v1);

            //! A public static member function.
            /*!
                y0(x)の適合点xfにおける関数値および微分値の型を変換する
                \param y y0(x)の適合点xfにおける関数値および微分値の型を変換する（std::array）
                \return y0(x)のx（原点に近い点）における関数値および微分値（Eigen::VectorXd）
            */
            static Eigen::VectorXd score(shootfunc::state_type const & y);

            //! A public static member function.
            /*!
                d[y0(x)]/dxの関数値および微分値を求める
                \param y y0(x)の関数値および微分値（入力）
                \param dydx d[y0(x)]/dxの関数値および微分値（出力）
                \param x y0(x)のx
            */
            static void rhs(shootfunc::state_type const & y, shootfunc::state_type & dydx, double const x);

            // #endregion メンバ関数

            // #region メンバ定数

            //! A public member variable (constant expression).
            /*!
                微分を差分化するための適当な刻み幅のベクトル
            */
            static auto constexpr DELV = 1.0E-14;

            //! A public member variable (constant expression).
            /*!
                y'0(0)の初期値
            */
            static auto constexpr V1 = -1.588076779;

            // #endregion 

            // #region 禁止されたコンストラクタ・メンバ関数

            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            shootfunc() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト（未使用）
            */
            shootfunc(shootfunc const & dummy) = delete;

            //! A public member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト
                \return コピー元のオブジェクト
            */
            shootfunc & operator=(shootfunc const & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }
}

#endif  // _SHOOTFUNC_H_
