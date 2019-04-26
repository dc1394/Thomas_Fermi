/*! \file readinputfile.cpp
    \brief インプットファイルの読み込みを行うクラスの実装
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

#include "readinputfile.h"
#include <iostream>                     // for std::cerr
#include <stdexcept>                    // for std::runtime_error
#include <utility>                      // for std::move
#include <boost/assert.hpp>             // for BPOOST_ASSERT
#include <boost/cast.hpp>               // for boost::numeric_cast
#include <boost/range/algorithm.hpp>    // for boost::find, boost::transform
#include <boost/tokenizer.hpp>          // for boost::tokenizer

namespace thomasfermi {
    // #region staticメンバ変数

    ci_string const ReadInputFile::CHEMICAL_NUMBER = "chemical.number";

    // #endregion staticメンバ変数

    // #region コンストラクタ

    ReadInputFile::ReadInputFile(std::pair<std::string, bool> const & arg) :
        PData([this] { return std::cref(pdata_); }, nullptr), 
        ifs_(std::get<0>(arg).c_str()),
        lineindex_(1),
        pdata_(std::make_shared<Data>())
    {
        pdata_->useomp_ = std::get<1>(arg);
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数
        
    void ReadInputFile::readFile()
    {
        if (!ifs_.is_open())
            throw std::runtime_error("インプットファイルが開けませんでした");

        auto const errorendfunc = []() {
            throw std::runtime_error("インプットファイルが異常です");
        };

        // 原子番号を読み込む
        if (!readAtom()) {
            errorendfunc();
        }

        // グリッドの最小値を読み込む
        readValue("grid.xmin", XMIN_DEFAULT, pdata_->xmin_);

        // グリッドの最大値を読み込む
        readValue("grid.xmax", XMAX_DEFAULT, pdata_->xmax_);

        // グリッドのサイズを読み込む
        readValue("grid.num", GRID_NUM_DEFAULT, pdata_->grid_num_);

        // 許容誤差を読み込む
        readValue("eps", EPS_DEFAULT, pdata_->eps_);
                
        // マッチングポイントを読み込む
        if (!readMatchPoint()) {
            errorendfunc();
        }
        
        // Gauss-Legendre積分の分点を読み込む
        readValue("gauss.legendre.integ", GAUSS_LEGENDRE_INTEG_DEFAULT, pdata_->gauss_legendre_integ_);
        
        // エネルギーを求めるときのGauss-Legendre積分の分点を読み込む
        readValue("gauss.legendre.integ.norm", GAUSS_LEGENDRE_INTEG_NORM_DEFAULT, pdata_->gauss_legendre_integ_norm_);
        
        // Iterationの最大ループ回数を読み込む
        readValue("iteration.maxIter", ITERATION_MAXITER_DEFAULT, pdata_->iteration_maxiter_);

        // Iterationの一次混合の重みを読み込む
        if (!readIterationMixingWeight()) {
            errorendfunc();
        }
        
        // Iterationの収束判定条件の値を読み込む
        readValue("iteration.criterion", ITERATION_CRITERION_DEFAULT, pdata_->iteration_criterion_);
    }
    
    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void ReadInputFile::errorMessage(ci_string const & s) const
    {
        std::cerr << "インプットファイルに" << s.c_str() << "の行が見つかりませんでした" << std::endl;
    }

    void ReadInputFile::errorMessage(std::int32_t line, ci_string const & s1, ci_string const & s2) const
    {
        std::cerr << "インプットファイルの[" << s1.c_str() << "]の行が正しくありません" << std::endl;
        std::cerr << line << "行目, 未知のトークン:" << s2.c_str() << std::endl;
    }

    std::pair<std::int32_t, std::optional<ReadInputFile::strvec>> ReadInputFile::getToken(ci_string const & article)
    {
        std::array<char, BUFSIZE> buf;
        ifs_.getline(buf.data(), BUFSIZE);
        ci_string const line(buf.data());

        // もし一文字も読めなかったら
        if (!ifs_.gcount()) {
            errorMessage(article);
            return std::make_pair(-1, std::nullopt);
        }

        // 読み込んだ行が空、あるいはコメント行でないなら
        if (!line.empty() && (line[0] != '#')) {
            // トークン分割
            using boost_char_sep = boost::char_separator<char>;
            using boost_tokenizer = boost::tokenizer<boost_char_sep>;
 
            auto const strs = std::string(line.c_str());
            boost_tokenizer tok(strs, boost_char_sep(" \r\t"));
            
            auto const len = std::distance(tok.begin(), tok.end());
            if (!len) {
                return std::make_pair(1, std::nullopt);
            }

            strvec tokens(len);
            boost::transform(
                    tok,
                    tokens.begin(),
                    [](auto && str) { return ci_string(str.c_str()); });
            
            auto const itr(tokens.begin());

            if (*itr != article) {
                errorMessage(lineindex_, article, *itr);
                return std::make_pair(-1, std::nullopt);
            }

            return std::make_pair(0, std::make_optional<strvec>(std::move(tokens)));
        }
        else {
            return std::make_pair(1, std::nullopt);
        }
    }

    bool ReadInputFile::readAtom()
    {
        // 原子の種類を読み込む
        auto const chemsym(readData(ReadInputFile::CHEMICAL_NUMBER));
        if (!chemsym) {
            return false;
        }
        else {
            pdata_->Z_ = std::stod(std::string(chemsym->c_str()));
            return true;
        }
    }

    std::optional<ci_string> ReadInputFile::readData(ci_string const & article)
    {
        for (; true; lineindex_++) {
            auto const ret(getToken(article));

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens(*(std::get<1>(ret)));

                // 読み込んだトークンの数がもし2個以外だったら
                if (tokens.size() != 2 || tokens[1].empty()) {
                    std::cerr << "インプットファイル" << lineindex_ << "行目の、[" << article.c_str() << "]の行が正しくありません" << std::endl;
                    return std::nullopt;
                }

                ++lineindex_;
                return *(++tokens.begin());
            }
            case 1:
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい!");
                break;
            }
        }
    }

    std::optional<ci_string> ReadInputFile::readData(ci_string const & article, ci_string const & def)
    {
        // グリッドを読み込む
        for (; true; lineindex_++) {
            auto const ret(getToken(article));

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens(*(std::get<1>(ret)));
                auto itr(++tokens.begin());
                ++lineindex_;

                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    // デフォルト値を返す
                    return def;
                    break;

                case 2:
                    return *itr == "DEFAULT" ? def : *itr;
                    break;

                default:
                    {
                        auto val(*itr);

                        if (val == "DEFAULT" || val[0] == '#') {
                            // デフォルト値を返す
                            return def;
                        }
                        else if ((*(++itr))[0] != '#') {
                            errorMessage(lineindex_ - 1, article, *itr);
                            return std::nullopt;
                        }

                        return val;
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

    std::optional<ci_string> ReadInputFile::readDataAuto(ci_string const & article)
    {
        for (; true; lineindex_++) {
            auto const ret(getToken(article));

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens(*(std::get<1>(ret)));
                auto itr(++tokens.begin());
                ++lineindex_;

                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    return std::nullopt;
                    break;

                case 2:
                    return (*itr == "DEFAULT" || *itr == "AUTO") ?
                        std::make_optional<ci_string>() : std::make_optional<ci_string>(*itr);
                    break;

                default:
                    {
                        auto val = *itr;

                        if (val == "DEFAULT" || val == "AUTO" || val[0] == '#') {
                            // デフォルト値を返す
                            return std::make_optional<ci_string>();
                        } else if ((*(++itr))[0] != '#') {
                            errorMessage(lineindex_ - 1, article, *itr);

                            // エラー
                            return std::nullopt;
                        }

                        return std::make_optional<ci_string>(std::move(val));
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
    
    bool ReadInputFile::readIterationMixingWeight()
    {
        readValue("iteration.Mixing.Weight", ITERATION_MIXING_WEIGHT_DEFAULT, pdata_->iteration_mixing_weight_);
        if (pdata_->iteration_mixing_weight_ <= 0.0 || pdata_->iteration_mixing_weight_ > 1.0) {
            std::cerr << "インプットファイルの[iteration.Mixing.Weight]の行が正しくありません" << std::endl;
            return false;
        }
        return true;
    }

    bool ReadInputFile::readMatchPoint()
    {
        readValue("matching.point", MATCH_POINT_DEFAULT, pdata_->match_point_);
        if (pdata_->match_point_ <= pdata_->xmin_ || pdata_->match_point_ >= pdata_->xmax_) {
            std::cerr << "インプットファイルの[matching.point]の行が正しくありません" << std::endl;
            return false;
        }
        return true;
    }

    
    // #endregion privateメンバ関数
}

