﻿/*! \file gr_pulay.cpp
	\brief guaranteed-reduction Pulay法でyの合成を行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
    (but this is originally adapted by T.Ozaki for GR_Pulay.c from ADPACK)
	This software is released under the BSD-2 License.
*/

#include "gr_pulay.h"
#include <boost/cast.hpp>	// boost::numeric_cast
#include <Eigen/Core>		// Eigen::MatrixXd
#include <Eigen/Dense>

namespace thomasfermi {
    namespace mixing {
        // #region コンストラクタ

        GR_Pulay::GR_Pulay(std::shared_ptr<Data> const & pdata) :
            SimpleMixing(pdata),
            ryarray_(GR_Pulay::NUM_MIXING_PDM),
            yarray_(GR_Pulay::NUM_MIXING_PDM)
        {
        }

        // #region コンストラクタ

        // #region publicメンバ関数

        femall::FEM::dmklvector GR_Pulay::operator()(std::int32_t scfiter, femall::FEM::dvector const & x, femall::FEM::dmklvector const & y)
        {
			yarray_[0] = y;

			switch (scfiter) {
			case 1:
				yarray_[1] = Yold;
				return gr(1, x);

			default:
			{
				if (scfiter <= GR_Pulay::NUM_MIXING_PDM) {
					return gr(1, x);
				}
				else {
					return gr(scfiter - (GR_Pulay::NUM_MIXING_PDM - 2), x);
				}
			}
			break;
			}
        }

        // #endregion publicメンバ関数
        
        // #region privateメンバ関数

        femall::FEM::dmklvector GR_Pulay::getry()
        {
            return getry(yarray_[0], yarray_[1]);
        }

        femall::FEM::dmklvector GR_Pulay::getry(femall::FEM::dmklvector const & newy, femall::FEM::dmklvector const & oldy)
        {
            auto const size = newy.size();
            BOOST_ASSERT(size == oldy.size());

            femall::FEM::dmklvector ry(size);
            for (auto i = 0U; i < size; i++) {
                ry[i] = newy[i] - oldy[i];
            }

            return ry;
        }

		femall::FEM::dmklvector GR_Pulay::gr(std::int32_t scfiter, femall::FEM::dvector const & x)
		{
			auto const size = x.size();
			BOOST_ASSERT(size == Yold().size());

			switch (scfiter) {
			case 1:
			case 2:
				setyryarray();
				return yarray_[0];

			default:
			{
				if (scfiter & 1) {
					// Py
					py = yarray_[0];

					// Calc of Ry1
					ryarray_[1] = getry();

					return yarray_[0];
				}
				else {
					// Calc of RDM0
					ryarray_[0] = getry(yarray_[0], py);

					auto const max = ((scfiter - 2) >> 1) + 1;

					std::int32_t nummix, numslide;
					if (max < GR_Pulay::NUM_MIXING_PDM) {
						nummix = max;
						numslide = nummix + 1;
					}
					else{
						nummix = GR_Pulay::NUM_MIXING_PDM;
						numslide = nummix;
					}

					// alpha from RDM
					Eigen::MatrixXd a(nummix + 1, nummix + 1);

					ryarray_.resize(nummix + 1);
					for (auto scfi = 0; scfi <= nummix; scfi++) {
						for (auto scfj = scfi; scfj <= nummix; scfj++) {
							a(scfi, scfj) = 0.0;

							for (auto i = 0; i < size; i++){
								a(scfi, scfj) += ryarray_[scfi][i] * ryarray_[scfj][i];
							}

							a(scfj, scfi) = a(scfi, scfj);
						}
					}

					auto const av_dia = a(0, 0);
					for (auto scfi = 0; scfi <= nummix; scfi++) {
						for (auto scfj = 0; scfj <= nummix; scfj++){
							a(scfi, scfj) /= av_dia;
						}
					}

					Eigen::MatrixXd const ia(a.inverse());

					auto denominator = 0.0;
					for (auto scfi = 0; scfi <= nummix; scfi++) {
						for (auto scfj = scfi; scfj <= nummix; scfj++) {
							denominator += ia(scfi, scfj);
						}
					}

					std::vector<double> alden(nummix + 1);
					for (auto scfi = 0; scfi <= nummix; scfi++) {
						auto numerator = 0.0;
						for (auto scfj = 0; scfj <= nummix; scfj++) {
							numerator += ia(scfj, scfi);
						}
						alden[scfi] = numerator / denominator;
					}

					// Calculate an optimized residual rho
					std::vector<double> optry(size);

					auto optnorm_ry = 0.0;
					for (auto i = 0U; i < size; i++) {
						double dx;
						if (!i) {
							dx = x[0];
						}
						else {
							dx = x[i] - x[i - 1];
						}
						
						optry[i] = 0.0;
						for (auto pscfiter = 0; pscfiter <= nummix; pscfiter++) {
							optry[i] += alden[pscfiter] * ryarray_[pscfiter][i];
						}

						optnorm_ry += optry[i] * optry[i] * dx;
					}

					double coef_optrdm;
					if (1.0e-2 <= optnorm_ry) {
						coef_optrdm = 0.2;
					}
					else if (1.0e-4 <= optnorm_ry && optnorm_ry < 1.0e-2) {
						coef_optrdm = 0.3;
					}
					else if (1.0e-6 <= optnorm_ry && optnorm_ry < 1.0e-4) {
						coef_optrdm = 0.4;
					}
					else if (1.0e-8 <= optnorm_ry && optnorm_ry < 1.0e-6) {
						coef_optrdm = 0.5;
					}
					else {
						coef_optrdm = 1.0;
					}

					// Mixing of DM
					for (auto i = 0U; i < size; i++) {
						yarray_[0][i] = 0.0;
						for (auto pscf_iter = 0; pscf_iter <= nummix; pscf_iter++) {
							if (!pscf_iter) {
								yarray_[0][i] += alden[pscf_iter] * py[i];
							}
							else {
								yarray_[0][i] += alden[pscf_iter] * yarray_[pscf_iter][i];
							}
						}

						// Correction by the optimized rho
						yarray_[0][i] += coef_optrdm * optry[i];
						if (yarray_[0][i] < 0.0) {
							yarray_[0][i] = 1.0e-15;
						}
					}

					if (yarray_.size() >= boost::numeric_cast<std::size_t>(numslide)) {
						yarray_.resize(numslide + 1);
					}

					// Shift of rho
					for (auto pscf_iter = numslide; 0 < pscf_iter; pscf_iter--) {
						yarray_[pscf_iter] = yarray_[pscf_iter - 1];
					}

					if (ryarray_.size() >= boost::numeric_cast<std::size_t>(numslide)) {
						ryarray_.resize(numslide + 1);
					}

					// Shift of residual rho
					for (auto pscf_iter = numslide; 1 < pscf_iter; pscf_iter--){
						ryarray_[pscf_iter] = ryarray_[pscf_iter - 1];
					}
				}

				return yarray_[0];
			}
			break;
			}
		}

		void GR_Pulay::setyryarray()
        {
            ryarray_[2] = getry();

            yarray_[0] = SimpleMixing::operator()(yarray_[0], yarray_[1]);
            yarray_[2] = yarray_[1];
            yarray_[1] = yarray_[0];
        }

        // #endregion privateメンバ関数
    }
}

