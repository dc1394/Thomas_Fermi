﻿/*! \file soelement.h
    \brief 二次要素のクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/


#include "soelement.h"

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        SOElement::SOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk)
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
                (*plnods_)[i][0] = 2 * i;
                (*plnods_)[i][1] = 2 * i + 2;
                (*plnods_)[i][2] = 2 * i + 1;
            }
        }
        
        // #endregion コンストラクタ

        // #region publicメンバ関数

        FEM::resulttuple SOElement::createresult() const
        {
            return std::make_tuple(a0_, a1_, a2_, b_);
        }
        
        void SOElement::reset(dvector const & beta)
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
                    auto const lnodi = (*plnods_)[ielem][i];
                    auto const lnodj = (*plnods_)[ielem][j];

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

        FEM::dvector SOElement::getdndr(double r) const
        {
            dvector dndr(ntnoel_);
            dndr[0] = r - 0.5;
            dndr[1] = r + 0.5;
            dndr[2] = - 2.0 * r;

            return std::move(dndr);
        }

        FEM::dvector SOElement::getc(std::size_t ielem) const
        {
            auto const xl = coords_[(*plnods_)[ielem][1]] - coords_[(*plnods_)[ielem][0]];

            dvector c(ntnoel_);
            c[0] = gl_.qgauss(
                myfunctional::make_functional(
                    [this, ielem](double r)
                   { return - N1_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]] + N3_(r) * coords_[(*plnods_)[ielem][2]]); }),
                   -1.0,
                   1.0) * xl * 0.5;
            c[1] = gl_.qgauss(
                myfunctional::make_functional([this, ielem](double r)
                   { return - N2_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]] + N3_(r) * coords_[(*plnods_)[ielem][2]]); }),
                   -1.0,
                   1.0) * xl * 0.5;
            c[2] = gl_.qgauss(
                myfunctional::make_functional([this, ielem](double r)
                   { return - N3_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]] + N3_(r) * coords_[(*plnods_)[ielem][2]]); }),
                   -1.0,
                   1.0) * xl * 0.5;

            return c;
        }

        // #endregion privateメンバ関数
    }
}
