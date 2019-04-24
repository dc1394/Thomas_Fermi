/*! \file foelement.cpp
    \brief 一次要素のクラスの実装
    Copyright © 2015-2019 @dc1394 All Rights Reserved.
    (but this is originally adapted by 渡辺浩志 for stiff5.c from http://www.sml.k.u-tokyo.ac.jp/members/nabe/FEM/FEM.pdf )
    
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
#include "myfunctional/functional.h"

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        FOElement::FOElement(std::vector<double> && beta, std::vector<double> const & coords, std::size_t nint, bool useomp) :
            FEM(std::move(beta), coords, nint, useomp),
            func_([this](double x) { return pbeta_->operator()<Element::First>(x); }),
            fun1_([this](double r, double xl, std::size_t ielem)
            { return -N1_(r) * func_(N1_(r) * coords_[lnods_[ielem][0]] + N2_(r) * coords_[lnods_[ielem][1]]) * xl * 0.5; }),
            fun2_([this](double r, double xl, std::size_t ielem)
            { return -N2_(r) * func_(N1_(r) * coords_[lnods_[ielem][0]] + N2_(r) * coords_[lnods_[ielem][1]]) * xl * 0.5; })
        {
            auto const N1tmp = [](double r) { return 0.5 * (1.0 - r); };
            N1_ = std::cref(N1tmp);
            
            auto const N2tmp = [](double r) { return 0.5 * (1.0 + r); };
            N2_ = std::cref(N2tmp);

            ntnoel_ = 2;
            nelem_ = nnode_ - 1;
            
            initialize();
            
            for (auto i = 0U; i < nelem_; i++) {
                lnods_[i][0] = i;
                lnods_[i][1] = i + 1;
            }
        }

        // #endregion コンストラクタ

        // #region publicメンバ関数
        
        FEM::resulttuple FOElement::createresult() const
        {
            return std::forward_as_tuple(a0_, a1_, std::vector<double>(), b_);
        }

        void FOElement::reset(std::vector<double> const & beta)
        {
            FEM::reset(beta);
            func_ = [this](double x) { return pbeta_->operator()<Element::First>(x); };
        }

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        void FOElement::amerge(std::size_t ielem)
        {
            a0_[ielem] += astiff_[0][0];
            a0_[ielem + 1] += astiff_[1][1];
            a1_[ielem] = astiff_[0][1];
        }

        void FOElement::element(std::size_t ielem)
        {
            astiffclear();

            for (auto ir = 0U; ir < nint_; ir++) {
                auto const dndr(getdndr());
                
                FEM::element(dndr, ielem, ir);
            }
        }

        std::vector<double> FOElement::getc(std::size_t ielem) const
        {
            auto const xl = coords_[lnods_[ielem][1]] - coords_[lnods_[ielem][0]];
            std::vector<double> c(ntnoel_);

            c[0] = gl_.qgauss(
                myfunctional::make_functional([this, xl, ielem](double r){ return fun1_(r, xl, ielem); }),
                -1.0,
                1.0);

            c[1] = gl_.qgauss(
                myfunctional::make_functional([this, xl, ielem](double r){ return fun2_(r, xl, ielem); }),
                -1.0,
                1.0);

            return c;
        }

        std::vector<double> FOElement::getdndr() const
        {
            std::vector<double> dndr(ntnoel_);
            dndr[0] = - 0.5;
            dndr[1] = 0.5;

            return dndr;
        }

        // #endregion privateメンバ関数
    }
}

