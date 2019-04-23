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

#ifndef _LENEAREQUATIONS_H_
#define _LENEAREQUATIONS_H_

#pragma once

#include "element.h"
#include <tuple>    // for std::tuple
#include <vector>   // for std::vector

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
                デストラクタ
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
    }
}

#endif  // _LENEAREQUATIONS_H_

