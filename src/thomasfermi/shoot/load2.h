/*! \file load2.h
    \brief y(x)の初期関数y0(x)の、原点と端点における
           関数値とその微分値を求めるクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _LOAD2_H_
#define _LOAD2_H_

#pragma once

#include "shootfunc.h"
#include "../utility/deleter.h"
#include <memory>				// for std::unique_ptr

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

            //! A constructor.
            /*!
                唯一のコンストラクタ
            */
            load2();

            //! A destructor.
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
			std::unique_ptr<gsl_interp_accel, decltype(utility::gsl_interp_accel_deleter)> const acc_;

			//! A private member variable.
			/*!
				gsl_splineへのスマートポインタ
			*/
			std::unique_ptr<gsl_spline, decltype(utility::gsl_spline_deleter)> const spline_;

            // #endregion メンバ変数
		};
	}
}

#endif  // _LOAD2_H_
