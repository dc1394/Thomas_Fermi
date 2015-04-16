#include "Iteration.h"
/*! \file iteration.cpp
	\brief 微分方程式を反復法で解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/


#include "shoot/shootf.h"
#include <iostream>
#include <dvec.h>
#include <boost/cast.hpp>
#include <boost/format.hpp>
#include <boost/assert.hpp>
#include <boost/assign.hpp>
#include <boost/utility/in_place_factory.hpp>

namespace thomasfermi {
	namespace fem_all {
		// #region コンストラクタ・デストラクタ

		Iteration::Iteration(double alpha, double dx, std::size_t n, double eps, bool usesimd, bool usetbb, double x1, double x2, double xf) :
			alpha_(alpha),
			eps_(eps),
			ple_(boost::none),
			usesimd_(usesimd),
			usetbb_(usetbb)
		{
			using namespace thomasfermi;
			using namespace thomasfermi::shoot;

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
				l2.make_v2(x2));
			
			shootf::result_type xyvtuple(s(x1, x2, xf));
			y1_ = std::get<1>(xyvtuple)[0];
			y2_ = std::get<1>(xyvtuple).back();
			v1_ = std::get<2>(xyvtuple)[0];
			
			x_ = std::get<0>(xyvtuple);
			y_ = ybefore_ = FEM::dmklvector(std::get<1>(xyvtuple).begin(), std::get<1>(xyvtuple).end());

			pfem_.reset(new fem_all::FOElement(n, useSSEorAVX_, usecilk_, x_, Iteration::make_beta()));
			pfem_->stiff();

			i_bc_given_.reserve(N_BC_GIVEN);

			using namespace boost::assign;

			i_bc_given_ += 0, (pfem_->getnnode() - 1);
			v_bc_nonzero_.reserve(N_BC_GIVEN);
			v_bc_nonzero_ += y1_, y2_;

			ple_ = boost::in_place(pfem_->createresult());

			ple_->bound(N_BC_GIVEN, i_bc_given_, N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

			y_ = ple_->LEsolver();
		}

		FEM::dvector Iteration::make_beta() const
		{
			const std::size_t size = y_.size();
			BOOST_ASSERT(size == x_.size());
			FEM::dvector beta(size);

			for (std::size_t i = 0; i < size; i++)
				beta[i] = y_[i] * std::sqrt(y_[i] / x_[i]);

			return std::move(beta);
		}

		void Iteration::Iterationloop()
		{
			std::int32_t cnt = 0;
			double scferr = IterationTHRESHOLD, scferrbefore;
			do {
				ymix();
				pfem_->reset(Iteration::make_beta());
				pfem_->stiff2();

				ple_->reset(pfem_->getb());
				ple_->bound(N_BC_GIVEN, i_bc_given_, N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

				ybefore_ = y_;
				y_ = ple_->LEsolver();
				scferrbefore = scferr;
				scferr = IterationError();
				if (scferr > scferrbefore)
					alpha_ *= IterationREDUCTION;

				cnt++;
				std::cout << "反復回数: " << cnt << "回, IterationError: " << boost::format("%.15f\n") % scferr;
			} while (scferr > TOL_);

			pbeta_ = pfem_->getpbeta();
		}

		void Iteration::ymix()
		{
			const std::size_t size = y_.size();
			BOOST_ASSERT(size == ybefore_.size());
			
			for (std::size_t i = 0; i < size; i++)
				y_[i] = ybefore_[i] + alpha_ * (y_[i] - ybefore_[i]);

		}

		double Iteration::IterationError() const
		{
			const std::size_t size = y_.size();
			BOOST_ASSERT(size == ybefore_.size());

			double sum = 0.0;
			for (std::size_t i = 0; i < size; i++)
				sum += sqr(y_[i] - ybefore_[i]);

			return std::sqrt(sum);
		}
	}
}
