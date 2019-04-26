/*! \file getcomlineoption.h
    \brief コマンドラインオプションの解析を行うクラスの宣言
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

#ifndef _GETCOMLINEOPTION_H_
#define _GETCOMLINEOPTION_H_

#pragma once

#include <cstdint>  // for std::int32_t
#include <string>   // for std::string
#include <utility>  // for std::pair

namespace thomasfermi {
    //! A class.
    /*!
        コマンドラインオプションを解析するクラス
    */
    class GetComLineOption final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A default constructor.
        /*!
            デフォルトコンストラクタ
        */
        GetComLineOption() = default;

        //! A default destructor.
        /*!
            デフォルトデストラクタ
        */
        ~GetComLineOption() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            コマンドラインオプションを解析する
            \param argc コマンドライン引数の数
            \param argv コマンドライン引数
            \return 解析に成功したら0、失敗したら-1
        */
        std::int32_t getopt(int argc, char * const argv[]);

        //! A public member function (constant).
        /*!
            インプットファイル名とCilk Plusを使用するかどうかを、std::pairで返す
            \return インプットファイル名とCilk Plusを使用するかどうかのstd::pair
        */
        std::pair<std::string, bool> getpairdata() const;

        // #endregion メンバ関数

    private:
        // #region メンバ変数

        //!  A private static member variable (constant).
        /*!
            デフォルトのインプットファイル名
        */
        static std::string const DEFINPNAME;

        //!  A private member variable.
        /*!
            デフォルトのインプットファイル名
        */
        std::string inpname_;

        //!  A private member variable.
        /*!
            OpenMPを使用するかどうか
        */
        bool useomp_ = false;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        GetComLineOption(GetComLineOption const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        GetComLineOption & operator=(GetComLineOption const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif  // _GETCOMLINEOPTION_H_
