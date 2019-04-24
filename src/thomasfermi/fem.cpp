/*! \file fem.h
    \brief 有限要素法のクラスの実装
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

#include "fem.h"
#include <cstdint>							// for std::uint32_t
#include <utility>							// for std::move
#include <boost/range/algorithm/fill.hpp>	// for boost::fill
#include <boost/assert.hpp>					// for BOOST_ASSERT
#include <omp.h>

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        FEM::FEM(std::vector<double> && beta, std::vector<double> const & coords, std::size_t nint, bool useomp) :
            B([this] { return std::cref(b_); }, nullptr),
            Nnode([this] { return nnode_; }, nullptr),
            Ntnoel([this] { return ntnoel_; }, nullptr),
            PBeta([this] { return std::cref(pbeta_); }, nullptr),
            nnode_(coords.size()),
            a0_(nnode_, 0.0),
            a1_(nnode_ - 1, 0.0),
            b_(nnode_, 0.0),
            beta_(std::move(beta)),
            coords_(coords),
            gl_(nint),
            nint_(nint),
            pbeta_(std::make_shared<Beta>(coords_, beta_)),
            useomp_(useomp)
        {
            BOOST_ASSERT(coords_.size() == beta_.size());
        }
        
        // #endregion コンストラクタ

        // #region publicメンバ関数

        void FEM::reset(std::vector<double> const & beta)
        {
            pbeta_.reset();
            pbeta_ = std::make_shared<Beta>(coords_, beta);

			boost::fill(b_, 0.0);
        }

        void FEM::stiff()
        {
            for (auto ielem = 0U; ielem < nelem_; ielem++) {
                element(ielem);
            }

            if (useomp_) {
                auto const nelem = static_cast<std::int32_t>(nelem_);
#pragma omp parallel for
                for (auto ielem = 0; ielem < nelem; ielem++) {
                    amerge(ielem);
                    createb(ielem);
                }
            }
            else {
                for (auto ielem = 0U; ielem < nelem_; ielem++) {
                    amerge(ielem);
                    createb(ielem);
                }
            }
        }

        void FEM::stiff2()
        {
            if (useomp_) {
                auto const nelem = static_cast<std::int32_t>(nelem_);
#pragma omp parallel for
                for (auto ielem = 0; ielem < nelem; ielem++) {
                    createb(ielem);
                }
            }
            else {
                for (auto ielem = 0U; ielem < nelem_; ielem++) {
                    createb(ielem);
                }
            }
        }

        // #endregion publicメンバ関数

        // #region protectedメンバ関数

        void FEM::astiffclear()
        {
            for (auto i = 0U; i < ntnoel_; i++) {
                for (auto j = 0U; j < ntnoel_; j++) {
                    astiff_[i][j] = 0.0;
                }
            }
        }

        void FEM::element(std::vector<double> const & dndr, std::size_t ielem, std::size_t ir)
        {
            auto ajacob = 0.0;

            for (auto i = 0U; i < ntnoel_; i++) {
                ajacob += dndr[i] * coords_[lnods_[ielem][i]];
            }

            auto const detjac = ajacob;
            auto const ajainv = 1.0 / ajacob;

            std::vector<double> dndx(ntnoel_);
            for (auto i = 0U; i < ntnoel_; i++) {
                dndx[i] = dndr[i] * ajainv;
            }

            auto const detwei = detjac * gl_.W()[ir];
            for (auto i = 0U; i < ntnoel_; i++) {
                for (auto j = 0U; j < ntnoel_; j++) {
                    astiff_[i][j] += detwei * dndx[i] * dndx[j];
                }
            }
        }

        void FEM::initialize()
        {
            lnods_.resize(boost::extents[nelem_][ntnoel_]);
        }

        // #endregion protectedメンバ関数

        // #region privateメンバ関数

        void FEM::createb(std::size_t ielem)
        {
            auto const c(getc(ielem));
            for (auto i = 0U; i < ntnoel_; i++) {
                b_[lnods_[ielem][i]] += c[i];
            }
        }

        // #endregion privateメンバ関数
    }
}

