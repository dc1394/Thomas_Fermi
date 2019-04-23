/*! \file checkpoint.h
    \brief 時間計測のためのクラスの宣言
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

#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

#pragma once

#include "fastarenaobject.h"
#include <array>				// for std::array			
#include <chrono>               // for std::chrono               
#include <cstdint>              // for std::int32_t, std::int64_t
#include <memory>               // for std::unique_ptr
#include <utility>              // for std::pair

namespace checkpoint {
    //! A class.
    /*!
        時間計測のためのクラス
    */
	class CheckPoint final {
        // #region クラスの前方宣言

        //! A structure.
        /*!
            チェックポイントの情報を格納する構造体
        */
	    struct Timestamp {
            //! A public member variable.
            /*!
                行数
            */
		    std::int32_t line;

            //! A public member variable.
            /*!
                チェックポイントの名称
            */
            char const * action;

            //! A public member variable.
            /*!
                チェックポイントの時間
            */
            std::chrono::high_resolution_clock::time_point realtime;
	    };
                
        //! A struct.
        /*!
            チェックポイントの情報の配列を格納する構造体
        */
	    struct CheckPointFastImpl {
            // #region コンストラクタ・デストラクタ

            //! A default constructor.
            /*!
                デフォルトコンストラクタかつ唯一のコンストラクタ
            */
            CheckPointFastImpl() : cur(0) {}

            //! A default destructor.
            /*!
                デフォルトデストラクタ
            */
            ~CheckPointFastImpl() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region メンバ変数

            //! A public static member variable (constant).
            /*!
                チェックポイントの数
            */
            static auto constexpr N = 30U;

            //! A public member variable.
            /*!
                現在の場所
            */
		    std::int32_t cur;
		
            //! A public member variable.
            /*!
                チェックポイントの情報の配列
            */
            std::array<CheckPoint::Timestamp, N> points;

            // #endregion メンバ変数
	    };

        // #endregion クラスの前方宣言

        // #region クラス内クラスの宣言と実装

        template <typename T>
        struct fastpimpl_deleter {
            void operator()(T * p) const {
                FastArenaObject<sizeof(CheckPointFastImpl)>::
                    operator delete(reinterpret_cast<void *>(p));
            }
        };

        // #endregion クラス内クラスの宣言と実装

    public:
        // #region コンストラクタ・デストラクタ

        //! A default constructor.
        /*!
            デフォルトコンストラクタかつ唯一のコンストラクタ
        */
        CheckPoint();

        //! A default destructor.
        /*!
            デフォルトデストラクタ
        */
        ~CheckPoint() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            チェックポイントを設定する
            \param line 行数
            \param action チェックポイントの名称
        */
        void checkpoint(char const * action, std::int32_t line);

        //! A public member function.
        /*!
            直前のチェックポイントから計測した、経過時間を表示する
        */
        void checkpoint_print() const;

        //! A public member function.
        /*!
            最初のチェックポイントから最後のチェックポイント
            までの経過時間を表示する
        */
        void totalpassageoftime() const;

        // #endregion メンバ関数 

    private:
        // #region メンバ変数

        //! A private member variable (constant).
        /*!
            チェックポイントの状態へのスマートポインタ
        */
        std::unique_ptr< CheckPointFastImpl, fastpimpl_deleter<CheckPointFastImpl> > const cfp;

        // #endregion メンバ変数

    public:
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト
        */
		CheckPoint(CheckPoint const & dummy) = delete;

        //! operator=() (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト
            \return コピー元のオブジェクト
        */
		CheckPoint & operator=(CheckPoint const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

    // #region 非メンバ関数

    //! A function.
    /*!
        自分自身のプロセスのメモリ使用量を計測する    
    */
    void usedmem();

    // #endregion 非メンバ関数
}

#endif  // _CHECKPOINT_H_

