/*! \file ci_string.h
    \brief 大文字小文字を区別しない文字列クラスの宣言
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

#ifndef _CI_STRING_H_
#define _CI_STRING_H_

#pragma once

#include <string>   // for std::string
#include <cctype>   // for std::toupper

namespace thomasfermi {
    //! A class.
    /*!
        大文字小文字を区別しない文字列クラス
    */
    struct ci_char_traits : public std::char_traits<char> {
        // #region staticメンバ関数

        //! A public static member function.
        /*!
            2個のバッファー (大文字と小文字を区別しない) を比較する
            \param s1 一つ目の文字列
            \param s2 二つ目の文字列
            \param 文字列の長さ
            \return std::strcmpの戻り値と同様
        */
        static int compare(char const * s1, char const * s2, std::size_t n);

        //! A public static member function.
        /*!
            引数で与えられた二つの文字が等しいかどうかを判別する
            \param lhs 一つ目の文字
            \param rhs 二つ目の文字
            \return 二つの文字が等しいかどうか
        */
        static bool eq(char lhs, char rhs)
        {
            return std::toupper(lhs) == std::toupper(rhs);
        }
        
        //! A public static member function.
        /*!
            引数で与えられた二つの文字の大小関係を判別する
            \param lhs 一つ目の文字
            \param rhs 二つ目の文字
            \return 二つの文字の大小関係
        */
        static bool lt(char lhs, char rhs)
        {
            return std::toupper(lhs) < std::toupper(rhs);
        }

        // #endregion staticメンバ関数
    };

    // #region 型エイリアス

    using ci_string = std::basic_string < char, ci_char_traits >;

    // #endregion 型エイリアス

    // #region 非メンバ関数

#ifndef _MSC_VER
    //! A function.
    /*!
        2個のバッファー (大文字と小文字を区別しない) を比較する
        \param s1 一つ目の文字列
        \param s2 二つ目の文字列
        \param 文字列の長さ
        \return std::strcmpの戻り値と同様
    */
    int memIcmp(char const * s1, char const * s2, std::size_t n);
#endif

    // #endregion 非メンバ関数
}

#endif  // _CISTRING_H_

