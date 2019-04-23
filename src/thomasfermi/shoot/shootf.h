/*! \file shootf.h
    \brief 狙い撃ち法により、y(x)を求めるクラスの宣言
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

#ifndef _SHOOTF_H_
#define _SHOOTF_H_

#pragma once

#include "load2.h"
#include "../myfunctional/functional.h"
#include <functional>					// for std::function
#include <utility>						// for std::pair
#include <vector>						// for std::vector

namespace thomasfermi {
    namespace shoot {
        //! A class.
        /*!
            狙い撃ち法により、y(x)を求めるクラス
        */
        class shootf final {
            // #region 型エイリアス

			using loadfunctype = std::function<shootfunc::state_type(double, double)>;

        public:
			using result_type = std::pair<std::vector<double>, std::vector<double>>;

        private:
			using scorefunctype = std::function<Eigen::VectorXd(shootfunc::state_type const &)>;

            // #endregion 型エイリアス

        public:
            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
                \param delv1 原点に近いxにおけるyの微分値の変化
                \param delv2 無限遠点に近いxにおけるyの微分値の変化
                \param dx xのメッシュの間隔
                \param eps 許容誤差
                \param load1 原点に近いxにおけるyの値とその微分値を求める関数オブジェクト
                \param l2 原点と無限遠点に近いxにおけるyの値とその微分値を求める関数オブジェクト
                \param score 適合点で合致するべきyの値とその微分値を求める関数オブジェクト
                \param v1 原点に近いxにおけるyの微分値
                \param v2 無限遠点に近いxにおけるyの微分値
            */
            shootf(double delv1,
                   double delv2,
                   double dx,
                   double eps,
                   loadfunctype const & load1,
                   load2 const & l2,
                   scorefunctype const & score,
                   double v1,
                   double v2);

            //! A default destructor.
            /*!
				デフォルトデストラクタ
            */
			~shootf() = default;

            // #endregion コンストラクタ・デストラクタ

        private:
            // #region privateメンバ関数

            //! A private member function (const).
            /*!
                最終的な結果を生成する
                \param res1 原点に近い点から適合点までの結果
                \param res2 無限遠点に近い点から適合点までの結果
                \param x1 原点に近いxの値
                \param xf 適合点に近いxの値
                \return xのメッシュとそれに対応したyの値のtuple
            */
			shootf::result_type createResult(std::vector<double> const & res1, std::vector<double> const & res2, double x1, double xf) const;

            // #endregion privateメンバ関数

            // #region publicメンバ関数
            
        public:
            //! A public member function (const).
            /*!
                結果を生成する
				\param usecilk cilkを使うかどうか
                \param x1 原点に近いxの値
                \param x2 無限遠点に近いxの値
                \param xf 適合点のxの値
                \return xのメッシュとそれに対応したyの値のtuple
            */
			shootf::result_type operator()(bool usecilk, double x1, double x2, double xf);

            // #endregion publicメンバ関数

        private:
            // #region メンバ変数

            //! A private member variable (constant expression).
            /*!
                許容誤差
            */
            static double constexpr EPS = 1.0E-14;

            //! A private member variable (constant).
            /*!
                原点に近いxにおけるyの微分値の増分
            */
            double const delv1_;

            //! A private member variable (constant).
            /*!
                無限遠点のxにおけるyの微分値の増分
            */
            double const delv2_;

            //! A private member variable (constant).
            /*!
                xのメッシュの間隔
            */
            double const dx_;

            //! A private member variable (constant).
            /*!
                許容誤差
            */
            double const eps_;

			//! A private member variable (constant).
			/*!
				原点に近いxにおけるyの値とその微分値を求める関数オブジェクト
			*/
			loadfunctype const load1_;

			//! A private member variable (constant).
			/*!
				無限遠点に近いxにおけるyの値とその微分値を求める関数オブジェクト
			*/
			loadfunctype const load2_;

			//! A private member variable (constant).
			/*!
				適合点で合致するべきyの値とその微分値を求める関数オブジェクト
			*/
			scorefunctype const score_;

            //! A private member variable.
            /*!
                原点に近いxにおけるyの微分値
            */
            double v1_;

            //! A private member variable.
            /*!
                無限遠点に近いxにおけるyの微分値
            */
            double v2_;

            // #region 禁止されたコンストラクタ・メンバ関数

        public:
            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            shootf() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト
            */
            shootf(shootf const & dummy) = delete;

            //! A public member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト
                \return コピー元のオブジェクト
            */
            shootf & operator=(shootf const & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
	}
}

#endif  // _SHOOTF_H_

