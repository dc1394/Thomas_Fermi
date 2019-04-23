/*! \file iteration.cpp
    \brief 微分方程式を反復法で解くクラスの実装
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

#include "foelement.h"
#include "iteration.h"
#include "readinputfile.h"
#include "shoot/shootf.h"
#include <iostream>                             // for std::cout
#include <stdexcept>                            // for std::runtime_error
#include <boost/assert.hpp>                     // for BOOST_ASSERT
#include <boost/format.hpp>                     // for boost::format

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ・デストラクタ

        Iteration::Iteration(std::pair<std::string, bool> const & arg) :
            PData([this] { return std::cref(pdata_); }, nullptr)
        {
            using namespace thomasfermi;
            using namespace thomasfermi::shoot;

            // インプットファイルの読み込み
            ReadInputFile rif(arg);         
            rif.readFile();
            pdata_ = rif.PData;

            // 混合法オブジェクトの生成
            pmix_ = std::make_unique<mixing::SimpleMixing>(pdata_);

            // メッシュの間隔を求める
            auto const dx = pdata_->xmax_ / static_cast<double>(pdata_->grid_num_);

            load2 l2;
            shootf s(
                shootfunc::DELV,
                shootfunc::DELV,
                dx,
                shootfunc::DELV * 0.1,
                shootfunc::load1,
                l2,
                shootfunc::score,
                shootfunc::V1,
                l2.make_v2(pdata_->xmax_));
            
            auto const usecilk = std::get<1>(arg);
            auto const [xtmp, ytmp] = s(usecilk, pdata_->xmin_, pdata_->xmax_, pdata_->match_point_);

            x_ = xtmp;
            y_ = std::vector<double>(ytmp.cbegin(), ytmp.cend());
            y1_ = ytmp[0];
            y2_ = ytmp.back();
            
            pmix_->Yold = y_;

            pfem_.reset(new femall::FOElement(make_beta(), x_, pdata_->gauss_legendre_integ_, usecilk));
            pfem_->stiff();

            i_bc_given_.reserve(Iteration::N_BC_GIVEN);

            i_bc_given_ = { 0, pfem_->Nnode - 1 };
            v_bc_nonzero_.reserve(Iteration::N_BC_GIVEN);
            v_bc_nonzero_ = { y1_, y2_ };

            ple_.emplace(pfem_->createresult());

            ple_->bound<Element::First>(Iteration::N_BC_GIVEN, i_bc_given_, Iteration::N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

            y_ = ple_->LEsolver<Element::First>();
        }

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数
        
        void Iteration::Iterationloop()
        {
            for (auto i = 1U; i < pdata_->iteration_maxiter_; i++) {
                pfem_->reset(Iteration::make_beta());
                pfem_->stiff2();

                ple_->reset(pfem_->B);
                ple_->bound<Element::First>(Iteration::N_BC_GIVEN, i_bc_given_, Iteration::N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

                pmix_->Yold = y_;

                ymix(ple_->LEsolver<Element::First>());

                auto const normrd = GetNormRD();

                std::cout << "反復回数: " << i << "回, NormRD: " << boost::format("%.15f\n") % normrd;
                if (normrd < pdata_->iteration_criterion_) {
                    pbeta_ = pfem_->PBeta;
                    return;
                }
            }
            
            throw std::runtime_error("収束しませんでした。");
        }

        Iteration::result_type Iteration::makeresult()
        {
            auto const y_prime_0 = (y_[1] - y_[0]) / (x_[1] - x_[0]);
            return std::forward_as_tuple(std::move(pbeta_), std::move(x_), y_prime_0);
        }

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        double Iteration::GetNormRD() const
        {
            auto const size = y_.size();

            auto const & yold(pmix_->Yold());

            BOOST_ASSERT(size == yold.size());

            auto sum = 0.0;
            for (auto i = 0U; i < size; i++) {
                sum += sqr(y_[i] - yold[i]);
            }

            return std::sqrt(sum);
        }

        std::vector<double> Iteration::make_beta() const
        {
            auto const size = y_.size();
            BOOST_ASSERT(size == x_.size());
            std::vector<double> beta(size);

            for (auto i = 0U; i < size; i++) {
                beta[i] = y_[i] * std::sqrt(y_[i] / x_[i]);
            }

            return beta;
        }

        void Iteration::ymix(std::vector<double> const & y)
        {
            y_ = (*pmix_)(y);
        }
        
        // #endregion privateメンバ関数
    }
}

