/*! \file fem.h
	\brief 有限要素法のクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	(but this is originally adapted by 渡辺浩志 for stiff5.c from http://www.sml.k.u-tokyo.ac.jp/members/nabe/FEM/FEM.pdf )
	This software is released under the BSD-2 License.
*/

#define _ITERATOR_DEBUG_LEVEL 0

#include "fem.h"
#include <cstdint>						// for std::uint32_t
#include <utility>						// for std::move
#include <boost/assert.hpp>				// for BOOST_ASSERT
#include <tbb/parallel_for.h>           // for tbb::parallel_for
#include <tbb/partitioner.h>            // for tbb::auto_partitioner
#include <tbb/task_scheduler_init.h>    // for tbb::task_scheduler_init

namespace thomasfermi {
	namespace fem_all {
		// #region コンストラクタ

		FEM::FEM(dvector && beta, dvector const & coords, std::size_t nint, bool usesimd, bool usetbb) :
			B([this] { return b_; }, nullptr),
			Nnode([this] { return nnode_; }, nullptr),
			Ntnoel([this] { return ntnoel_; }, nullptr),
			PBeta([this] { return pbeta_; }, nullptr),
			a1_(nnode_, 0.0),
			a2_(nnode_ - 1, 0.0),
			b_(nnode_, 0.0),
			beta_(std::move(beta)),
			coords_(coords),
			func_([this](double x) { return (*pbeta_)(x); }),
			gl_(nint),
			nint_(nint),
			nnode_(coords.size()),
			pbeta_(std::make_shared<Beta>(coords_, beta_)),
			usesimd_(usesimd),
			usetbb_(usetbb)
		{
			BOOST_ASSERT(coords_.size() == beta_.size());
		}
		
		// #endregion コンストラクタ

		// #region publicメンバ関数

		std::tuple<FEM::dmklvector, FEM::dmklvector, FEM::dmklvector> FEM::createresult() const
		{
			return std::make_tuple(a1_, a2_, b_);
		}

		void FEM::reset(const dvector & beta)
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

			if (usetbb_) {
				tbb::task_scheduler_init init;

				tbb::parallel_for(
					std::uint32_t(0),
					nelem_,
					std::uint32_t(1),
					[this](std::uint32_t ielem)
				{
					amerge(ielem);

					dvector const c(getc(ielem));
					for (auto i = 0U; i < ntnoel_; i++)
					{
						b_[lnods_[i][ielem]] += c[i];
					}
				},
					tbb::auto_partitioner());
			}
			else {
				for (auto ielem = 0U; ielem < nelem_; ielem++) {
					amerge(ielem);

					dvector const c(getc(ielem));
					for (auto i = 0U; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				}
			}
		}

		void FEM::stiff2()
		{
			if (usetbb_) {
				tbb::task_scheduler_init init;

				tbb::parallel_for(
					std::uint32_t(0),
					nelem_,
					std::uint32_t(1),
					[this](std::uint32_t ielem)
				{
					dvector const c(getc(ielem));
					for (auto i = 0U; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				},
					tbb::auto_partitioner());
			}
			else {
				for (auto ielem = 0U; ielem < nelem_; ielem++) {
					dvector const c(getc(ielem));
					for (auto i = 0U; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				}
			}
		}

		// #endregion publicメンバ関数

		// #region protectedメンバ関数

		void FEM::initialize()
		{
			astiff_.resize(boost::extents[ntnoel_][ntnoel_]);
			lnods_.resize(boost::extents[ntnoel_][nelem_]);
		}

		// #endregion protectedメンバ関数

		// #region privateメンバ関数

		void FEM::amerge(std::size_t ielem)
		{
			a1_[ielem] += astiff_[0][0];
			a1_[ielem + 1] += astiff_[1][1];
			a2_[ielem] = astiff_[0][1];
		}

		void FEM::element(std::size_t ielem)
		{
			for (auto i = 0U; i < ntnoel_; i++)
				for (auto j = 0U; j < ntnoel_; j++)
					astiff_[i][j] = 0.0;

			for (auto ir = 0U; ir < nint_; ir++) {
				const dvector dndr(getdndr());
				double ajacob = 0.0;

				for (auto i = 0U; i < ntnoel_; i++) {
					ajacob += dndr[i] * coords_[lnods_[i][ielem]];
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
						astiff_[i][j] += detwei * dndx[i] * dndx[j];
					}
				}
			}
		}

		// #endregion privateメンバ関数
	}
}
