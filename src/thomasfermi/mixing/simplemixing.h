/*! \file simplemixing.h
	\brief 一次混合法でyの合成を行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
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
                \param pdata データオブジェクト
            */
            SimpleMixing(std::shared_ptr<Data> const & pdata);

            //! A destructor.
            /*!
                デストラクタ
            */
            virtual ~SimpleMixing() = default;

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

        protected:
            //!  A private member variable.
            /*!
                前回のループのyの値
            */
            femall::FEM::dmklvector yold_;

            // #endregion privateメンバ変数

            // #region protectedメンバ変数

        protected:
            //!  A protected member variable.
            /*!
                データオブジェクト
            */
            std::shared_ptr<Data> pdata_;

            // #endregion protectedメンバ変数

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
