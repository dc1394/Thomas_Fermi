/*! \file simplemixing.h
	\brief 一次混合を行うクラスの宣言

Copyright ©  2015 @dc1394 All Rights Reserved.
This software is released under the BSD-2 License.
*/

#ifndef _SIMPLEMIXING_H_
#define _SIMPLEMIXING_H_

#pragma once

#include "../data.h"
#include "../fem.h"
#include <vector>

namespace thomasfermi {
	namespace mixing {
        class SimpleMixing {
            // #region コンストラクタ・デストラクタ

        public:
            //! A constructor.
            /*!
                唯一のコンストラクタ
            */
            SimpleMixing(std::shared_ptr<Data> const & pdata);

            //! A destructor.
            /*!
                デストラクタ
            */
            virtual ~SimpleMixing() noexcept;

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ関数

            //! A public member function.
            /*!
                一次混合法によって、yの合成を行う関数
            */
            femall::FEM::dmklvector operator()(femall::FEM::dmklvector const & y);

            // #endregion メンバ関数

            // #region メンバ変数

            //!  A private member variable.
            /*!
                データオブジェクト
            */
            std::shared_ptr<Data> pdata_;

            //! A private member variable.
            /*!
                前回のループのyの値の可変長配列
            */
            femall::FEM::dmklvector ybefore_;

            // #endregion メンバ変数

            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            SimpleMixing() = delete;

            //! A private copy constructor (deleted).
            /*!
            コピーコンストラクタ（禁止）
            */
            SimpleMixing(const SimpleMixing &) = delete;

            //! A private member function (deleted).
            /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
            */
            SimpleMixing & operator=(const SimpleMixing &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }
}

#endif	// _SIMPLEMIXING_H_
