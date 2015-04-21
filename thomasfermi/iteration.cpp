﻿/*! \file iteration.cpp
	\brief 微分方程式を反復法で解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "iteration.h"
#include "shoot/shootf.h"
#include <iostream>								// for std::cout
#include <boost/cast.hpp>						// for boost::numeric_cast
#include <boost/format.hpp>						// for boost::format
#include <boost/assert.hpp>						// for BOOST_ASSERT
#include <boost/utility/in_place_factory.hpp>	// for boost::in_place

namespace thomasfermi {
	namespace femall {
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
			
			shootf::result_type xytuple(s(x1, x2, xf));
			y1_ = std::get<1>(xytuple)[0];
			y2_ = std::get<1>(xytuple).back();
			
			x_ = std::get<0>(xytuple);
			y_ = ybefore_ = FEM::dmklvector(std::get<1>(xytuple).begin(), std::get<1>(xytuple).end());

			pfem_.reset(new femall::FOElement(make_beta(), x_, n, usesimd_, usetbb_));
			pfem_->stiff();

			i_bc_given_.reserve(Iteration::N_BC_GIVEN);

			i_bc_given_ = { 0, pfem_->Nnode - 1 };
			v_bc_nonzero_.reserve(Iteration::N_BC_GIVEN);
			v_bc_nonzero_ = { y1_, y2_ };

			ple_ = boost::in_place(pfem_->createresult());

			ple_->bound(Iteration::N_BC_GIVEN, i_bc_given_, Iteration::N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

			y_ = ple_->LEsolver();
		}

		Iteration::~Iteration()
		{
			ple_ = boost::none;
		}

		// #endregion コンストラクタ・デストラクタ

		// #region publicメンバ関数
		
		void Iteration::Iterationloop()
		{
			auto cnt = 0;
			auto scferr = Iteration::ITERATION_THRESHOLD;
			double scferrbefore;
			do {
				ymix();
				pfem_->reset(Iteration::make_beta());
				pfem_->stiff2();

				ple_->reset(pfem_->B);
				ple_->bound(Iteration::N_BC_GIVEN, i_bc_given_, Iteration::N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

				ybefore_ = y_;
				y_ = ple_->LEsolver();
				scferrbefore = scferr;
				scferr = IterationError();

				if (scferr > scferrbefore) {
					alpha_ *= Iteration::ITERATION_REDUCTION;
				}

				cnt++;
				std::cout << "反復回数: " << cnt << "回, IterationError: " << boost::format("%.15f\n") % scferr;
			} while (scferr > eps_);

			pbeta_ = pfem_->PBeta;
		}

		Iteration::result_type Iteration::makeresult()
		{
			return std::make_pair(std::move(pbeta_), std::move(x_));
		}

		// #endregion publicメンバ関数

		// #region privateメンバ関数

		double Iteration::IterationError() const
		{
			auto const size = y_.size();
			BOOST_ASSERT(size == ybefore_.size());

			auto sum = 0.0;
			for (auto i = 0U; i < size; i++) {
				sum += sqr(y_[i] - ybefore_[i]);
			}

			return std::sqrt(sum);
		}

		FEM::dvector Iteration::make_beta() const
		{
			auto const size = y_.size();
			BOOST_ASSERT(size == x_.size());
			FEM::dvector beta(size);

			for (auto i = 0U; i < size; i++) {
				beta[i] = y_[i] * std::sqrt(y_[i] / x_[i]);
			}

			return std::move(beta);
		}

		void Iteration::ymix()
		{
			auto const size = y_.size();
			BOOST_ASSERT(size == ybefore_.size());

			for (auto i = 0U; i < size; i++) {
				y_[i] = ybefore_[i] + alpha_ * (y_[i] - ybefore_[i]);
			}
		}
		
		// #endregion privateメンバ関数
	}
}