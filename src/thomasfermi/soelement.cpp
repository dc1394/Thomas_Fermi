/*! \file soelement.h
    \brief 二次要素のクラスの実装
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


#include "soelement.h"

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        SOElement::SOElement(std::vector<double> && beta, std::vector<double> const & coords, std::size_t nint, bool usecilk)
            :   FEM(std::move(beta), coords, nint, usecilk),
                a2_(nnode_ - 2, 0.0),
                func_([this](double x) { return pbeta_->operator() < Element::Second > (x); })
        {
            auto const N1tmp = [](double r) { return -0.5 * r * (1.0 - r); };
            N1_ = std::cref(N1tmp);

            auto const N2tmp = [](double r) { return 0.5 * r * (1.0 + r); };
            N2_ = std::cref(N2tmp);
            
            auto const N3tmp = [](double r) { return (1.0 - r * r); };
            N3_ = std::cref(N3tmp);

            ntnoel_ = 3;
            nelem_ = (nnode_ - 1) >> 1;
            
            initialize();

            for (auto i = 0U; i < nelem_; i++) {
                lnods_[i][0] = 2 * i;
                lnods_[i][1] = 2 * i + 2;
                lnods_[i][2] = 2 * i + 1;
            }
        }
        
        // #endregion コンストラクタ

        // #region publicメンバ関数

        FEM::resulttuple SOElement::createresult() const
        {
            return std::forward_as_tuple(a0_, a1_, a2_, b_);
        }
        
        void SOElement::reset(std::vector<double> const & beta)
        {
            FEM::reset(beta);
            func_ = [this](double x) { return pbeta_->operator()< Element::Second >(x); };
        }

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        void SOElement::amerge(std::size_t ielem)
        {
            for (auto i = 0UL; i < ntnoel_; i++) {
                for (auto j = 0UL; j < ntnoel_; j++) {
                    auto const lnodi = lnods_[ielem][i];
                    auto const lnodj = lnods_[ielem][j];

                    if (lnodj == lnodi) {
                        a0_[lnodj] += astiff_[i][j];
                    }
                    else if (lnodj == lnodi - 1 && lnodj < nnode_ - 1) {
                        a1_[lnodj] += astiff_[i][j];
                    }
                    else if (lnodj == lnodi - 2 && lnodj < nnode_ - 2) {
                        a2_[lnodj] += astiff_[i][j];
                    }
                }
            }
        }
        
        void SOElement::element(std::size_t ielem)
        {
            astiffclear();

            for (auto ir = 0U; ir < nint_; ir++) {
                auto const dndr(getdndr(gl_.X()[ir]));

                FEM::element(dndr, ielem, ir);
            }
        }

        std::vector<double> SOElement::getdndr(double r) const
        {
            std::vector<double> dndr(ntnoel_);
            dndr[0] = r - 0.5;
            dndr[1] = r + 0.5;
            dndr[2] = - 2.0 * r;

            return dndr;
        }

        std::vector<double> SOElement::getc(std::size_t ielem) const
        {
            auto const xl = coords_[lnods_[ielem][1]] - coords_[lnods_[ielem][0]];

            std::vector<double> c(ntnoel_);
            c[0] = gl_.qgauss(
                myfunctional::make_functional(
                    [this, ielem](double r)
                   { return - N1_(r) * func_(N1_(r) * coords_[lnods_[ielem][0]] + N2_(r) * coords_[lnods_[ielem][1]] + N3_(r) * coords_[lnods_[ielem][2]]); }),
                   -1.0,
                   1.0) * xl * 0.5;
            c[1] = gl_.qgauss(
                myfunctional::make_functional([this, ielem](double r)
                   { return - N2_(r) * func_(N1_(r) * coords_[lnods_[ielem][0]] + N2_(r) * coords_[lnods_[ielem][1]] + N3_(r) * coords_[lnods_[ielem][2]]); }),
                   -1.0,
                   1.0) * xl * 0.5;
            c[2] = gl_.qgauss(
                myfunctional::make_functional([this, ielem](double r)
                   { return - N3_(r) * func_(N1_(r) * coords_[lnods_[ielem][0]] + N2_(r) * coords_[lnods_[ielem][1]] + N3_(r) * coords_[lnods_[ielem][2]]); }),
                   -1.0,
                   1.0) * xl * 0.5;

            return c;
        }

        // #endregion privateメンバ関数
    }
}

