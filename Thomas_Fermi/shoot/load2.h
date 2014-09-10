/*! \file load2.h
    \brief y(x)の初期関数y0(x)の、原点と端点における
           関数値とその微分値を求めるクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#ifndef _LOAD2_H_
#define _LOAD2_H_

#pragma once

#include "shootfunc.h"
#include "../Spline.h"
#include <boost/optional.hpp>

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
            */
            ~load2() { psplint = boost::none; }

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ関数
            
            //! A public static member function.
            /*!
                y0(x)のx（遠点）における微分値の近似値を求める
                \param x xの値
                \return dy0(x) / dxの値
            */
            static double dy0(double x);

            //! A public static member function.
            /*!
                y0(x)のx（遠点）における近似値を求める
                \param x xの値
                \return y0(x)の値
            */
            static double y0(double x);

            //! A public member function (const).
            /*!
            y0(x)の端点における微分値を求める
            \param x2 端点
            \return y0(x)の端点における微分値
            */
            shootfunc::tmpary make_v2(double x2) const;

            //! A public member function (const).
            /*!
                y0(x)の端点における値を求め、引数で与えた微分値との
                std::arrayを作って返す
                \param v2 y0(x)の端点における微分値
                \param x2 y0(x)の端点における値
                \return 端点における値とその微分値が収納されたstd::array
            */
            shootfunc::state_type operator()(shootfunc::tmpary const & v2, double x2) const;

            //! A public member variable (constant expression).
            /*!
                無限遠点における関数値とその微分値を求めるときに使う定数λ
            */
            static double constexpr Lambda = 3.886;

            //! A public member variable (constant expression).
            /*!
                補間と近似表式を切り替える時の閾値
            */
            static double constexpr Threshold = 65.0;

            //! A public member variable (constant expression).
            /*!
                無限遠点における関数値とその微分値を求めるときに使う定数k
            */
            static double constexpr K = 0.190785707092222;

            //! A public member variable (constant expression).
            /*!
                補間のために使う動的配列のサイズ
            */
			static std::size_t constexpr Xysize = 150;
                        
            //! A public member variable (constant expression).
            /*!
                Spline補間のために使う
            */
            boost::optional<Spline> psplint;

		};
	}
}

#endif  // _LOAD2_H_
