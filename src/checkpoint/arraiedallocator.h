﻿/*! \file arraiedallocator.h
    \brief 固定サイズのメモリを確保するアロケータークラス
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

#ifndef _ARRAYIEDALLOCATOR_H_
#define _ARRAYIEDALLOCATOR_H_

#pragma once

#include <cstdint>                  // for std::uint32_t
#include <boost/static_assert.hpp>  // for BOOST_STATIC_ASSERT

namespace checkpoint {
    //! A template class.
    /*!
        固定サイズのメモリを確保するアロケータークラス
        \param TTypeSize 収納する型のサイズ
        \param TnumArray 収納する要素の数
    */
    template <std::size_t TTypeSize, std::size_t TNumArray>
    class ArraiedAllocator final
    {
        // サイズは絶対０より大きくなくちゃダメ
        BOOST_STATIC_ASSERT(TNumArray > 0);

        // #region クラス内クラスの宣言と実装

        //! A structure.
        /*!
        */
        struct Item {
            union {
                char value_[TTypeSize];
                struct Item * next_;
            };
        };

        // #endregion クラス内クラスの宣言と実装 

    public:
        // #region コンストラクタ・デストラクタ

        //! A default constructor.
        /*!
            デフォルトコンストラクタかつ唯一のコンストラクタ
        */
        ArraiedAllocator();

        //! A default destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ArraiedAllocator() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public static member function.
        /*!
            メモリを確保してそのアドレスを返す
            \return 確保されたメモリのアドレス
        */
        static void * Alloc() {
            Item * ret = first_;
            first_ = ret->next_;
            return reinterpret_cast<void *>(ret);
        }

        //! A public static member function.
        /*!
            確保されたメモリを解放する
            \param item 解放するメモリのアドレス
        */
        static void Free(void * item) {
            Item * rev = reinterpret_cast<Item *>(item);
            rev->next_ = first_;
            first_ = rev;
        }

        //! A public static member function.
        /*!
            アロケーターを返す
            \return アロケーター
        */
        static ArraiedAllocator& GetAllocator() { return allocator_; }

        //! A public static member function.
        /*!
            格納できる最大の要素数を返す
            \return 収納できる最大の要素数
        */
        static std::size_t Max() { return MAX_SIZE; }

    private:
        // #region メンバ変数

        //! A private static member variable (constant expression).
        /*!
            格納できる最大の要素数
        */
        static const std::size_t MAX_SIZE = TNumArray;

        //! A private static member variable.
        /*!
            アロケーター
        */
        static ArraiedAllocator allocator_;

        //! A private static member variable.
        /*!
            最初の要素へのポインタ
        */
        static Item * first_;
        
        //! A private static member variable.
        /*!
            要素の配列
        */
        static Item items_[MAX_SIZE];

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        ArraiedAllocator(ArraiedAllocator const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト（未使用）
        */
        ArraiedAllocator & operator=(ArraiedAllocator const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    template <std::size_t TTypeSize, std::size_t TNumArray>
    inline ArraiedAllocator<TTypeSize, TNumArray>::ArraiedAllocator() {
        first_ = &items_[0];
        for (std::uint32_t i = 0; i < TNumArray; i++) {
            items_[i].next_ = &items_[i + 1];
        }
        items_[TNumArray - 1].next_ = nullptr;
    }

    template <std::size_t TTypeSize, std::size_t TNumArray>
    ArraiedAllocator<TTypeSize, TNumArray> ArraiedAllocator<TTypeSize, TNumArray>::allocator_;

    template <std::size_t TTypeSize, std::size_t TNumArray>
    typename ArraiedAllocator<TTypeSize, TNumArray>::Item* ArraiedAllocator<TTypeSize, TNumArray>::first_;

    template <std::size_t TTypeSize, std::size_t TNumArray>
    typename ArraiedAllocator<TTypeSize,TNumArray>::Item
        ArraiedAllocator<TTypeSize, TNumArray>::items_[ArraiedAllocator<TTypeSize,TNumArray>::MAX_SIZE];
}

#endif // _ARRAYIEDALLOCATOR_H_
