/*! \file fastarenaobject.h
    \brief 指定された型の指定された要素数のメモリを確保するクラス
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

#ifndef _FASTARENAOBJECT_H_
#define _FASTARENAOBJECT_H_

#pragma once

#include "arraiedallocator.h"

namespace checkpoint {
    //! A template class.
    /*!
        指定された型の指定された要素数のメモリを確保するクラス
        \param TTypeSize 収納する型のサイズ
        \param TnumArray 収納する要素の数
    */
	template <size_t TTypeSize, size_t TNumArray = 1>
	struct FastArenaObject final
	{
		// サイズは絶対０より大きくなくちゃダメ
		BOOST_STATIC_ASSERT(TNumArray > 0);

        // #region メンバ関数

        //! A public member function.
        /*!
            operator newの宣言と実装
            \param dummy 未使用
        */
		static void * operator new(std::size_t dummy) {
			return ArraiedAllocator<TTypeSize, TNumArray>::GetAllocator().Alloc();
		}

        //! A public member function.
        /*!
            operator deleteの宣言と実装
            \param p 解放するメモリの先頭アドレス
        */
		static void operator delete(void * p) {
			ArraiedAllocator<TTypeSize, TNumArray>::GetAllocator().Free(p);
		}

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A default constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        FastArenaObject() = delete;

        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        FastArenaObject(FastArenaObject const & dummy) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト（未使用）
        */
        FastArenaObject & operator=(FastArenaObject const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};
}

#endif  // _FASTARENAOBJECT_H_

