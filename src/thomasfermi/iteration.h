/*! \file iteration.h
    \brief 微分方程式を反復法で解くクラスの宣言
    Copyright © 2015-2019 @dc1394 All Rights Reserved.

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

#ifndef _ITERATION_H_
#define _ITERATION_H_

#pragma once

#include "data.h"
#include "fem.h"
#include "linearequations.h"
#include "mixing/simplemixing.h"
#include "shoot/shootfunc.h"
#include "utility/property.h"
#include <optional>                 // for std::nullopt, std::optional

namespace thomasfermi {
    namespace femall {
        class Iteration final {
            // #region 型エイリアス

        public:
            using result_type = std::tuple<std::shared_ptr<Beta>, std::vector<double>, double>;

            // #endregion 型エイリアス

            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
                唯一のコンストラクタ
                \param arg インプットファイル名とCilk Plusを使用するかどうかのstd::pair
            */
            explicit Iteration(std::pair<std::string, bool> const & arg);

            //! A default destructor.
            /*!
                デフォルトデストラクタ
            */
            ~Iteration()
            {
                ple_ = std::nullopt;
            }

            // #region コンストラクタ・デストラクタ

            // #region publicメンバ関数

            //! A public member function.
            /*!
                反復関数
            */
            void Iterationloop();

            //! A public member function.
            /*!
                結果を返す関数
                \return 結果
            */
            result_type makeresult();
            
            // #endregion publicメンバ関数

            // #region privateメンバ関数

            //! A private member function (const).
            /*!
                反復の誤差を返す
                \return 反復の誤差
            */
            double GetNormRD() const;
            
            //! A private member function (const).
            /*!
                βを生成する関数
                \return β
            */
            std::vector<double> make_beta() const;

            //! A private member function.
            /*!
                yを合成する
                \param scfiter SCFの回数
                \param y 新しいy
            */
            void ymix(std::vector<double> const & y);
                        
            // #endregion privateメンバ関数

            // #region プロパティ

        public:
            //! A property.
            /*!
                読み込んだデータへのプロパティ
            */
            utility::Property<std::shared_ptr<Data> const &> const PData;

            // #endregion プロパティ

            // #region メンバ変数

        private:
            //! A private member variable (constant expression).
            /*!
            */
            static auto constexpr N_BC_GIVEN = 2U;
            
            //! A private member variable.
            /*!
            */
            std::vector<std::size_t> i_bc_given_;

            //! A private member variable.
            /*!
                βオブジェクト
            */
            std::shared_ptr<Beta> pbeta_;
            
            //!  A private member variable.
            /*!
                データオブジェクト
            */
            std::shared_ptr<Data> pdata_;

            //! A private member variable.
            /*!
                有限要素法オブジェクト
            */
            std::unique_ptr<FEM> pfem_;

			//! A private member variable.
			/*!
				連立一次方程式のソルバーオブジェクト
			*/
            std::optional<Linear_equations> ple_;

            //! A private member variable.
            /*!
                yの混合法
            */
            std::unique_ptr<mixing::SimpleMixing> pmix_;

            //! A private member variable.
            /*!
            */
            std::vector<double> v_bc_nonzero_;

            //! A private member variable.
            /*!
                xのメッシュの可変長配列
            */
            std::vector<double> x_;
            
            //! A private member variable.
            /*!
                yの値の可変長配列
            */
            std::vector<double> y_;

            //! A private member variable.
            /*!
                yの境界条件（原点に近い方）
            */
            double y1_;

            //! A private member variable.
            /*!
                yの境界条件（原点から遠い方）
            */
            double y2_;

            // #endregion メンバ変数
                        
            // #region 禁止されたコンストラクタ・メンバ関数

        public:
            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            Iteration() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト（未使用）
            */
            Iteration(Iteration const & dummy) = delete;

            //! A public member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            Iteration & operator=(Iteration const & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }

    template <typename T>
    //! A template function（非メンバ関数）.
    /*!
        対象の値を二乗する
        \param x 対象の値
        \return 対象の値の二乗
    */
    constexpr inline T sqr(T x) noexcept
    {
        return x * x;
    }
}

#endif  // _ITERATION_H_

