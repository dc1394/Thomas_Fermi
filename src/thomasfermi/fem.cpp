/*! \file fem.h
	\brief 有限要素法のクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	(but this is originally adapted by 渡辺浩志 for stiff5.c from http://www.sml.k.u-tokyo.ac.jp/members/nabe/FEM/FEM.pdf )
	This software is released under the BSD-2 License.
*/

#include "fem.h"
#include <cstdint>				        // for std::uint32_t
#include <utility>				        // for std::move
#include <boost/assert.hpp>		        // for BOOST_ASSERT
#include <cilk/cilk.h>					// for cilik_for

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ

		FEM::FEM(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk) :
			B([this] { return b_; }, nullptr),
			Nnode([this] { return nnode_; }, nullptr),
			Ntnoel([this] { return ntnoel_; }, nullptr),
			PBeta([this] { return pbeta_; }, nullptr),
			nnode_(coords.size()),
			a1_(nnode_, 0.0),
			a2_(nnode_ - 1, 0.0),
			b_(nnode_, 0.0),
			beta_(std::move(beta)),
			coords_(coords),
			func_([this](double x) { return (*pbeta_)(x); }),
			gl_(nint),
			nint_(nint),
			pbeta_(std::make_shared<Beta>(coords_, beta_)),
			usecilk_(usecilk)
		{
			BOOST_ASSERT(coords_.size() == beta_.size());
		}
		
		// #endregion コンストラクタ

		// #region publicメンバ関数

		std::tuple<FEM::dmklvector, FEM::dmklvector, FEM::dmklvector> FEM::createresult() const
		{
			return std::make_tuple(a1_, a2_, b_);
		}

		void FEM::reset(dvector const & beta)
		{
			pbeta_.reset();
			pbeta_ = std::make_shared<Beta>(coords_, beta);
			func_ = [this](double x) { return (*pbeta_)(x); };

			for (double & v : b_) {
				v = 0.0;
			}
		}

		void FEM::stiff()
		{
			for (auto ielem = 0U; ielem < nelem_; ielem++)
				element(ielem);

			if (usecilk_) {
				cilk_for (auto ielem = 0U; ielem < nelem_; ielem++) {
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
			if (usecilk_) {
				cilk_for (auto ielem = 0U; ielem < nelem_; ielem++) {
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

		void FEM::initialize()
		{
			pastiff_.reset(new dmatrix(boost::extents[ntnoel_][ntnoel_]));
			plnods_.reset(new boost::multi_array<std::size_t, 2>(boost::extents[ntnoel_][nelem_]));
		}

		// #endregion protectedメンバ関数

		// #region privateメンバ関数

		void FEM::amerge(std::size_t ielem)
		{
			a1_[ielem] += (*pastiff_)[0][0];
			a1_[ielem + 1] += (*pastiff_)[1][1];
			a2_[ielem] = (*pastiff_)[0][1];
		}

        void FEM::createb(std::size_t ielem)
        {
            auto const c(getc(ielem));
            for (auto i = 0U; i < ntnoel_; i++) {
                b_[(*plnods_)[i][ielem]] += c[i];
            }
        }

		void FEM::element(std::size_t ielem)
		{
			for (auto i = 0U; i < ntnoel_; i++)
				for (auto j = 0U; j < ntnoel_; j++)
					(*pastiff_)[i][j] = 0.0;

			for (auto ir = 0U; ir < nint_; ir++) {
				auto const dndr(getdndr());
				auto ajacob = 0.0;

				for (auto i = 0U; i < ntnoel_; i++) {
					ajacob += dndr[i] * coords_[(*plnods_)[i][ielem]];
				}

				auto const detjac = ajacob;
				auto const ajainv = 1.0 / ajacob;
				
				dvector dndx(ntnoel_);
				for (auto i = 0U; i < ntnoel_; i++) {
					dndx[i] = dndr[i] * ajainv;
				}

				auto const detwei = detjac * gl_.W()[ir];
				for (auto i = 0U; i < ntnoel_; i++) {
					for (auto j = 0U; j < ntnoel_; j++) {
						(*pastiff_)[i][j] += detwei * dndx[i] * dndx[j];
					}
				}
			}
		}

		// #endregion privateメンバ関数
	}
}
