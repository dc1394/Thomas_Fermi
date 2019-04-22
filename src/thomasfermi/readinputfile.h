/*! \file readinputfile.h
    \brief インプットファイルの読み込みを行うクラスの宣言
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

#ifndef _READINPUTFILE_H_
#define _READINPUTFILE_H_

#pragma once

#include "ci_string.h"
#include "data.h"
#include "utility/property.h"
#include <fstream>                  // for std::ifstream
#include <memory>                   // for std::shared_ptr
#include <optional>                 // for std::optional
#include <vector>                   // for std::vector
#include <boost/lexical_cast.hpp>   // for boost::lexical_cast

namespace thomasfermi {
    //! A class.
    /*!
        インプットファイルを読み込み、Dataクラスのオブジェクトに格納するクラス
    */
    class ReadInputFile final {
        // #region 型エイリアス

        using strvec = std::vector < ci_string > ;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param arg インプットファイル名と、Cilk Plusを使用するかどうかのstd::pair
        */
        explicit ReadInputFile(std::pair<std::string, bool> const & arg);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ReadInputFile() = default;

        // #region コンストラクタ・デストラクタ
        
        // #region メンバ関数

        //! A public member function.
        /*!
            ファイルを読み込む
        */
        void readFile();

    private:
        //! A private member function (const).
        /*!
            エラーを表示する
            \param s エラーのトークン
            \return 読み込みが成功したかどうか
        */
        void errorMessage(ci_string const & s) const;

        //! A private member function (const).
        /*!
            エラーを表示する
            \param line エラーのある行
            \param s1 エラーのトークン
            \param s2 エラーのトークン2
        */
        void errorMessage(std::int32_t line, ci_string const & s, ci_string const & s2) const;

        //! A private member function (const).
        /*!
            解析対象の文字列をトークンに分割する
            \param article 解析対象の文字列
            \return 関数が成功したかどうかと、トークンのstd::pair
        */
        std::pair< std::int32_t, std::optional<ReadInputFile::strvec> > getToken(ci_string const & article);

        //! A private member function.
        /*!
            原子に関するデータを読み込む
            \return 読み込みが成功したかどうか
        */
        bool readAtom();

        //! A private member function.
        /*!
            文字列を解析して、データとして読み込んで返す
            \param article 解析対象の文字列
            \return 読みこんだ文字列データ
        */
        std::optional<ci_string> readData(ci_string const & article);

        //! A private member function.
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \param def デフォルトの文字列
            \return 読みこんだ文字列
        */
        std::optional<ci_string> readData(ci_string const & article, ci_string const & def);

        template <typename T>
        //! A private member function (template function).
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \param default_value デフォルト値
            \return 読みこんだ文字列
        */
        std::optional<T> readData(ci_string const & article, T const & default_value);

        //! A private member function.
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \return 読みこんだ文字列（読み込みに失敗したならboost::none）
        */
        std::optional<ci_string> readDataAuto(ci_string const & article);

        //! A private member function.
        /*!
            マッチングポイントの値を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readMatchPoint();

        //! A private member function.
        /*!
            Iterationの一次混合の重みを読み込む
            \return 読み込みが成功したかどうか
        */
        bool readIterationMixingWeight();

        template <typename T>
        //! A private member function.
        /*!
            対象の要素の値をその行から読み込む
            \param article 要素名
            \param default_value デフォルトの値
            \param value 読み込んだ値
        */
        void readValue(ci_string const & article, T const & default_value, T & value);

        template <typename T>
        //! A private member function.
        /*!
            対象の要素の値をその行から読み込む
            \param article 要素名
            \param value 読み込んだ値
            \return 読み込みが成功したかどうか
        */
        bool readValueAuto(ci_string const & article, std::optional<T> & value);

        // #endregion メンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            読み込んだデータを返す
            \return 読み込んだデータ
        */
        utility::Property<std::shared_ptr<Data> const &> const PData;

        // #endregion プロパティ

        // #region メンバ変数

    private:
        //! A private member variable (constant expression).
        /*!
            バッファサイズ
        */
        static std::streamsize constexpr BUFSIZE = 1024U;

        //! A private member variable (constant).
        /*!
            「chemical.number」の文字列
        */
        static ci_string const CHEMICAL_NUMBER;

        //! A private member variable.
        /*!
            ファイル読み込み用のストリーム
        */
        std::ifstream ifs_;
        
        //! A private member variable.
        /*!
            現在の行数
        */
        std::size_t lineindex_;

        //! A private member variable.
        /*!
            インプットファイルから読み込んだデータ
        */
        std::shared_ptr<Data> pdata_;
            
        // #endregion メンバ変数
        
        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A default constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ReadInputFile() = delete;

        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        ReadInputFile(ReadInputFile const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ReadInputFile & operator=(ReadInputFile const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    template <typename T>
    std::optional<T> ReadInputFile::readData(ci_string const & article, T const & default_value)
    {
        for (; ; lineindex_++) {
            auto const ret = getToken(article);
            switch (ret.first)
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens = *(ret.second);
                auto itr(++tokens.begin());

                lineindex_++;
                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    // デフォルト値を返す
                    return std::make_optional<T>(default_value);
                    break;

                case 2:
                    if (*(++itr) == "DEFAULT") {
                        // デフォルト値を返す
                        return std::make_optional<T>(default_value);
                    }
                    else {
                        try {
                            return std::make_optional<T>(boost::lexical_cast<T>(itr->c_str()));
                        }
                        catch (boost::bad_lexical_cast const &) {
                            errorMessage(lineindex_ - 1, article, *itr);
                            return std::nullopt;
                        }
                    }

                default:
                {
                    auto val = *itr;

                    if (val == "DEFAULT" || val[0] == '#') {
                        return std::make_optional<T>(default_value);
                    }
                    else if ((*(++itr))[0] != '#') {
                        errorMessage(lineindex_ - 1, article, *itr);
                        return std::nullopt;
                    }

                    try {
                        return std::make_optional<T>(boost::lexical_cast<T>(val.c_str()));
                    }
                    catch (boost::bad_lexical_cast const &) {
                        errorMessage(lineindex_ - 1, article, val);
                        return std::nullopt;
                    }
                }
                }
            }
                break;

            case 1:
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい!");
                break;
            }
        }
    }

    template <typename T>
    void ReadInputFile::readValue(ci_string const & article, T const & default_value, T & value)
    {
        auto const p = readData(article, default_value);
        if (p) {
            value = *p;
        }
        else {
            throw std::runtime_error("インプットファイルが異常です");
        }
    }

    template <typename T>
    bool ReadInputFile::readValueAuto(ci_string const & article, std::optional<T> & value)
    {
        if (auto const val = readDataAuto(article)) {
            if (!val->empty()) {
                try {
                    std::size_t idx;
                    auto const v = std::stod(val->c_str(), &idx);
                    if (idx != val->length()) {
                        throw std::invalid_argument("");
                    }
                    value = std::make_optional<double>(v);
                }
                catch (std::invalid_argument const &) {
                    errorMessage(lineindex_ - 1, article, *val);
                    return false;
                }
            }
            else {
                value = std::nullopt;
            }
        }
        else {
            return false;
        }

        return true;
    }
}

#endif  // _READINPUTFILE_H_

