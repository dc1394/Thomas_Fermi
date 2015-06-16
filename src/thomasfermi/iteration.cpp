/*! \file iteration.cpp
	\brief 微分方程式を反復法で解くクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "iteration.h"
#include "readinputfile.h"
#include "shoot/shootf.h"
#include <iostream>								// for std::cout
#include <stdexcept>							// for std::runtime_error
#include <boost/cast.hpp>						// for boost::numeric_cast
#include <boost/format.hpp>						// for boost::format
#include <boost/assert.hpp>						// for BOOST_ASSERT
#include <boost/utility/in_place_factory.hpp>	// for boost::in_place

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ・デストラクタ

		Iteration::Iteration(std::pair<std::string, bool> const & arg) :
			PData([this] { return pdata_; }, nullptr)
		{
			using namespace thomasfermi;
			using namespace thomasfermi::shoot;

			ReadInputFile rif(arg);         // ファイルを読み込む
			rif.readFile();
			pdata_ = rif.PData;
            pmix_ = std::make_unique<mixing::SimpleMixing>(pdata_);

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
			auto const xytuple(s(usecilk, pdata_->xmin_, pdata_->xmax_, pdata_->match_point_));
			y1_ = std::get<1>(xytuple)[0];
			y2_ = std::get<1>(xytuple).back();
			
			x_ = std::get<0>(xytuple);
			y_ = FEM::dmklvector(std::get<1>(xytuple).begin(), std::get<1>(xytuple).end());
            pmix_->Ybefore = y_;

			pfem_.reset(new femall::FOElement(make_beta(), x_, pdata_->gauss_legendre_integ_, usecilk));
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
			auto normrd = Iteration::ITERATION_THRESHOLD;
			double normrdbefore;
			for (auto i = 1U; i < pdata_->iteration_maxiter_; i++) {
				pfem_->reset(Iteration::make_beta());
				pfem_->stiff2();

				ple_->reset(pfem_->B);
				ple_->bound(Iteration::N_BC_GIVEN, i_bc_given_, Iteration::N_BC_GIVEN, i_bc_given_, v_bc_nonzero_);

                pmix_->Ybefore = y_;
				//y_ = ;
                ymix(i, ple_->LEsolver());
				normrdbefore = normrd;
				normrd = GetNormRD();

				/*if (normrd > normrdbefore) {
					pdata_->iteration_mixing_weight_ *= Iteration::ITERATION_REDUCTION;
				}*/

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
			return std::make_pair(std::move(pbeta_), std::move(x_));
		}

		// #endregion publicメンバ関数

		// #region privateメンバ関数

		double Iteration::GetNormRD() const
		{
			auto const size = y_.size();
            BOOST_ASSERT(size == pmix_->Ybefore().size());

			auto sum = 0.0;
			for (auto i = 0U; i < size; i++) {
                sum += sqr(y_[i] - pmix_->Ybefore()[i]);
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

        void Iteration::ymix(std::int32_t scfiter, femall::FEM::dmklvector const & y)
		{
            y_ = (*pmix_)(y);
		}
		
		// #endregion privateメンバ関数
	}
}
