/*! \file simplemixing.h
	\brief 一次混合法でyの合成を行うクラスの宣言
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

#ifndef _SIMPLEMIXING_H_
#define _SIMPLEMIXING_H_

#pragma once

#include "../data.h"
#include "../fem.h"
#include <vector>

namespace thomasfermi {
	namespace mixing {
        class SimpleMixing final {
            // #region コンストラクタ・デストラクタ

        public:
            //! A constructor.
            /*!
                唯一のコンストラクタ
                \param pdata データオブジェクト
            */
            SimpleMixing(std::shared_ptr<Data> const & pdata);

            //! A destructor.
            /*!
                デフォルトデストラクタ
            */
            ~SimpleMixing() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ関数

            //! A public member function.
            /*!
                一次混合法によって、yの合成を行う関数
                \param newy 合成前のy
                \return 合成後のy
            */
            femall::FEM::dmklvector operator()(femall::FEM::dmklvector const & y);

            // #endregion メンバ関数

            // #region プロパティ

        public:
            //! A property.
            /*!
                前回のループのyの値の可変長配列へのプロパティ
            */
            utility::Property<femall::FEM::dmklvector const &> Yold;

            // #endregion プロパティ

            // #region privateメンバ変数

        private:
            //!  A private member variable.
            /*!
                前回のループのyの値
            */
            femall::FEM::dmklvector yold_;

            //!  A protected member variable.
            /*!
                データオブジェクト
            */
            std::shared_ptr<Data> pdata_;

            // #endregion privateメンバ変数

            // #region 禁止されたコンストラクタ・メンバ関数

        public:
            //! A default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            SimpleMixing() = delete;

            //! A copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param dummy コピー元のオブジェクト（未使用）
            */
            SimpleMixing(const SimpleMixing &) = delete;

            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param dummy コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            SimpleMixing & operator=(const SimpleMixing & dummy) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }
}

#endif	// _SIMPLEMIXING_H_

