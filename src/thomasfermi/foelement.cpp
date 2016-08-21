﻿/*! \file foelement.cpp
    \brief 一次要素のクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    (but this is originally adapted by 渡辺浩志 for stiff5.c from http://www.sml.k.u-tokyo.ac.jp/members/nabe/FEM/FEM.pdf )
    This software is released under the BSD 2-Clause License.
*/

#include "foelement.h"
#include "myfunctional/functional.h"

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        FOElement::FOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk) :
            FEM(std::move(beta), coords, nint, usecilk),
            func_([this](double x) { return pbeta_->operator() < Element::First > (x); }),
            fun1_([this](double r, double xl, std::size_t ielem)
            { return -N1_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]]) * xl * 0.5; }),
            fun2_([this](double r, double xl, std::size_t ielem)
            { return -N2_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]]) * xl * 0.5; })
        {
            auto const N1tmp = [](double r) { return 0.5 * (1.0 - r); };
            N1_ = std::cref(N1tmp);
            
            auto const N2tmp = [](double r) { return 0.5 * (1.0 + r); };
            N2_ = std::cref(N2tmp);

            ntnoel_ = 2;
            nelem_ = nnode_ - 1;
            
            initialize();
            
            for (auto i = 0U; i < nelem_; i++) {
                (*plnods_)[i][0] = i;
                (*plnods_)[i][1] = i + 1;
            }
        }

        // #endregion コンストラクタ

        // #region publicメンバ関数
        
        FEM::resulttuple FOElement::createresult() const
        {
            return std::forward_as_tuple(a0_, a1_, FEM::dmklvector(), b_);
        }

        void FOElement::reset(dvector const & beta)
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

        FEM::dvector FOElement::getc(std::size_t ielem) const
        {
            auto const xl = coords_[(*plnods_)[ielem][1]] - coords_[(*plnods_)[ielem][0]];
            dvector c(ntnoel_);

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

        FEM::dvector FOElement::getdndr() const
        {
            dvector dndr(ntnoel_);
            dndr[0] = - 0.5;
            dndr[1] = 0.5;

            return dndr;
        }

        // #endregion privateメンバ関数
    }
}
