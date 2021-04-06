/*! \file linearequations.h
    \brief 連立方程式を解くクラスの宣言
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

#ifndef _LINEAREQUATIONS_H_
#define _LINEAREQUATIONS_H_

#pragma once

#include "element.h"
#include <cmath>                // for std::abs
#include <stdexcept>            // for std::logic_error, std::invalid_argument
#include <tuple>                // for std::tuple
#include <vector>               // for std::vector
#include <boost/format.hpp>     // for boost::format
#include <mkl_lapack.h>         // for dptsv_, dpbsv_

namespace thomasfermi {
    namespace femall {
        //! A class.
        /*!
            連立方程式を解くクラス
        */
        class Linear_equations final {
            // #region 型エイリアス

            using sivector = std::vector<std::size_t>;

            // #endregion 型エイリアス

            // #region コンストラクタ・デストラクタ

        public:
            //! A constructor.
            /*!
                唯一のコンストラクタ
                \param res 連立方程式のベクトル
            */
            explicit Linear_equations(FEM::resulttuple const & res);

            //! A default destructor.
            /*!
                デフォルトデストラクタ
            */
            ~Linear_equations() = default;

            // #endregion コンストラクタ・デストラクタ 

            // #region publicメンバ関数

            template <Element E>
            //! A public member function.
            /*!
                境界条件を定める
                \param n_bc_given 既知量の数
                \param i_bc_given 既知量のインデックス
                \param n_bc_nonzero 非零の既知量の数
                \param i_bc_nonzero 非零の既知量のインデックス
                \param v_bc_nonzero 非零の既知量
            */
            void bound(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero); 
            
            //! A public member function.
            /*!
                ベクトルbを初期化する
                \param b 対象のベクトルb
            */
            void reset(std::vector<double> const & b);
            
            template <Element E>
            //! A public member function.
            /*!
                連立一次方程式の解を求める
                \return 連立一次方程式の解
            */
            std::vector<double> LEsolver();

            // #endregion publicメンバ関数

            // #region メンバ変数

        private:
            //! A private member variable.
            /*!
                ベクトルa0
            */
            std::vector<double> a0_;

            //! A private member variable (constant).
            /*!
                ベクトルa0の複製
            */
            std::vector<double> a0back_;

            //! A private member variable.
            /*!
                ベクトルa1
            */
            std::vector<double> a1_;

            //! A private member variable (constant).
            /*!
                ベクトルa1の複製
            */
            std::vector<double> a1back_;

            //! A private member variable.
            /*!
                ベクトルa2
            */
            std::vector<double> a2_;

            //! A private member variable.
            /*!
                ベクトルb
            */
            std::vector<double> b_;

            //! A private member variable (constant).
            /*!
                ベクトルの要素数
            */
            std::size_t const n_;

            // #endregion メンバ変数        
            
            // #region 禁止されたコンストラクタ・メンバ関数

        public:
            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            Linear_equations() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト（未使用）
            */
            Linear_equations(Linear_equations const & dummy) = delete;

            //! A public member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            Linear_equations & operator=(Linear_equations const & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };

        // #region templateメンバ関数の実装

        template <>
        inline void Linear_equations::bound<Element::First>(std::size_t n_bc_given, sivector const & i_bc_given, std::size_t n_bc_nonzero, sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
        {
            b_[i_bc_nonzero[0] + 1] -= v_bc_nonzero[0] * a1_[i_bc_nonzero[0]];
            b_[i_bc_nonzero[1] - 1] -= v_bc_nonzero[1] * a1_[i_bc_nonzero[1] - 1];

            for (auto i = 0U; i < n_bc_nonzero; i++) {
                b_[i_bc_nonzero[i]] = v_bc_nonzero[i];
            }

            for (auto i = 0U; i < n_bc_given; i++) {
                a0_[i_bc_given[i]] = 1.0;
                a1_[i_bc_given[i] - i] = 0.0;
            }
        }

        template <>
        inline void Linear_equations::bound<Element::Second>(std::size_t n_bc_given, Linear_equations::sivector const & i_bc_given, std::size_t n_bc_nonzero, Linear_equations::sivector const & i_bc_nonzero, std::vector<double> const & v_bc_nonzero)
        {
            // Dirichlet boundary condition
            Linear_equations::bound<Element::First>(n_bc_given, i_bc_given, n_bc_nonzero, i_bc_nonzero, v_bc_nonzero);

            b_[i_bc_nonzero[0] + 2] -= v_bc_nonzero[0] * a2_[i_bc_nonzero[0]];
            b_[i_bc_nonzero[1] - 2] -= v_bc_nonzero[1] * a2_[i_bc_nonzero[1] - 2];

            a2_[i_bc_given[0]] = 0.0;
            a2_[i_bc_given[1] - 2] = 0.0;
        }

        template <>
        inline std::vector<double> Linear_equations::LEsolver<Element::First>()
        {
            auto n = static_cast<std::int32_t>(n_);
            auto nrhs = 1;

            std::int32_t info;
            dptsv_(
                &n,
                &nrhs,
                a0_.data(),
                a1_.data(),
                b_.data(),
                &n,
                &info);

            if (info > 0) {
                throw std::logic_error("U is singular");
            }
            else if (info < 0) {
                auto const str = (boost::format("%d-th argument has illegal value") % std::abs(info)).str();

                throw std::invalid_argument(str);
            }

            return b_;
        }

        template <>
        inline std::vector<double> Linear_equations::LEsolver<Element::Second>()
        {
            // 上三角要素を使う場合
            auto uplo = 'U';

            // 線形方程式の数（行列Aの次数）
            auto n = static_cast<std::int32_t>(n_);

            // 係数行列の帯の中にある対角線より上の部分の個数
            auto kd = 2;

            // 行列{B}の列数。通常通り1
            auto nrhs = 1;

            // 配列ABの1次元目の大きさ（= KD + 1）
            auto nb = kd + 1;

            // 係数行列の帯の外を省略して詰め込んだ2次元配列
            // ピボッティングありのLU分解を行うために(KD + 1)× N必要
            std::vector<double> ab(nb * n);

            for (auto i = 0; i < n; i++) {
                for (auto j = i; j <= i + 2; j++) {
                    if (j == i) {
                        ab[(j) * nb + (kd + i - j)] = a0_[i];
                    }
                    else if (j == i + 1 && j < n - 1) {
                        ab[(j) * nb + (kd + i - j)] = a1_[i];
                    }
                    else if (j == i + 2 && j < n - 2) {
                        ab[(j) * nb + (kd + i - j)] = a2_[i];
                    }
                }
            }

            std::int32_t info;
            dpbsv_(
                &uplo,      // 上三角要素を使う場合                    
                &n,         // 線形方程式の数（行列Aの次数）
                &kd,        // 係数行列の帯の中にある対角線より上の部分の個数
                &nrhs,      // 行列{B}の列数。通常通り1
                ab.data(),  // 係数行列(input)，コレスキー分解の結果(output)
                &nb,        // 配列ABの1次元目の大きさ（=KD+1） 
                b_.data(),  // 方程式の右辺(input)，方程式の解(output)
                &n,         // 行列Bの1次元目の大きさ（=N）
                &info);

            if (info > 0) {
                throw std::logic_error("U is singular");
            }
            else if (info < 0) {
                auto const str = (boost::format("%d-th argument has illegal value") % std::abs(info)).str();

                throw std::invalid_argument(str);
            }

            return b_;
        }

        // #endregion templateメンバ関数の実装
    }
}

#endif  // _LINEAREQUATIONS_H_
