/*! \file load2.h
    \brief y(x)の初期関数y0(x)の、原点と端点における
           関数値とその微分値を求めるクラスの宣言
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

#ifndef _LOAD2_H_
#define _LOAD2_H_

#pragma once

#include "shootfunc.h"
#include <memory>				// for std::unique_ptr
#include <gsl/gsl_spline.h>     // for gsl_interp_accel, gsl_interp_accel_free, gsl_spline, gsl_spline_free

namespace thomasfermi {
	namespace shoot {
        //! A class.
        /*!
            y(x)の初期関数y0(x)の、原点と端点における
            関数値とその微分値を求めるクラス
        */
		class load2	final {
        public:
            // #region コンストラクタ

            //! A default constructor.
            /*!
                デフォルトコンストラクタかつ唯一のコンストラクタ
            */
            load2();

            //! A default destructor.
            /*!
				デフォルトデストラクタ
            */
			~load2() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ関数
            
            //! A public static member function.
            /*!
                y0(x)のx（遠点）における微分値の近似値を求める
                \param x xの値
                \return dy0(x) / dxの値
            */
            static double dy0(double x) noexcept;

            //! A public static member function.
            /*!
                y0(x)のx（遠点）における近似値を求める
                \param x xの値
                \return y0(x)の値
            */
            static double y0(double x) noexcept;

            //! A public member function (const).
            /*!
            y0(x)の端点における微分値を求める
            \param x2 端点
            \return y0(x)の端点における微分値
            */
            double make_v2(double x2) const noexcept;

            //! A public member function (const).
            /*!
                y0(x)の端点における値を求め、引数で与えた微分値との
                std::arrayを作って返す
                \param v2 y0(x)の端点における微分値
                \param x2 y0(x)の端点における値
                \return 端点における値とその微分値が収納されたstd::array
            */
            shootfunc::state_type operator()(double v2, double x2) const;

            // #endregion メンバ関数

            // #region メンバ変数

            //! A public member variable (constant expression).
            /*!
                無限遠点における関数値とその微分値を求めるときに使う定数λ
            */
            static auto constexpr LAMBDA = 3.886;

            //! A public member variable (constant expression).
            /*!
                補間と近似表式を切り替える時の閾値
            */
            static auto constexpr THRESHOLD = 65.0;

            //! A public member variable (constant expression).
            /*!
                無限遠点における関数値とその微分値を求めるときに使う定数k
            */
            static auto constexpr K = 0.190785707092222;

            //! A public member variable (constant expression).
            /*!
                補間のために使う動的配列のサイズ
            */
			static auto constexpr XYSIZE = 150;

        private:
			//! A private member variable.
			/*!
				gsl_interp_accelへのスマートポインタ
			*/
			std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> const acc_;

			//! A private member variable.
			/*!
				gsl_splineへのスマートポインタ
			*/
			std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)> const spline_;

            // #endregion メンバ変数
            
            // #region 禁止されたコンストラクタ・メンバ関数

        public:
            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト（未使用）
            */
			load2(load2 const & dummy) = delete;
			
            //! A public member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト
                \return コピー元のオブジェクト
            */
            load2 & operator=(load2 const & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif  // _LOAD2_H_

