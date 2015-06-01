/*! \file gr_pulay.cpp
	\brief guaranteed-reduction Pulay法でyの合成を行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
    (but this is originally adapted by T.Ozaki for GR_Pulay.c from ADPACK)
	This software is released under the BSD-2 License.
*/

#include "gr_pulay.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace thomasfermi {
    namespace mixing {
        // #region コンストラクタ

        GR_Pulay::GR_Pulay(std::shared_ptr<Data> const & pdata) :
            SimpleMixing(pdata),
            ryarray_(GR_Pulay::NUM_MIXING_PDM)
        {
        }

        // #region コンストラクタ

        // #region publicメンバ関数

        femall::FEM::dmklvector GR_Pulay::operator()(std::int32_t scfiter, femall::FEM::dvector const & x, femall::FEM::dmklvector const & y)
        {
            auto const size = y.size();
            BOOST_ASSERT(size == Ybefore.size());

            switch (scfiter) {
            case 0:
                yarray_[0] = y;
                yarray_[1] = Ybefore;

                setyryarray();
                return yarray_[0];
                break;

            case 1:
                setyryarray();
                return yarray_[0];

            default:
            {
                femall::FEM::dmklvector py;
                std::int32_t nummix, numslide;

                if (!(scfiter & 1)) {
                    // Py
                    py = yarray_[0];
                    
                    // Calc of Ry1
                    ryarray_[1] = getry();
                }
                else {
                    // Calc of RDM0
                    ryarray_[0] = getry(ryarray_[0], py);

                    auto const max = ((scfiter - 1) >> 1) + 1;
                    
                    if (max < GR_Pulay::NUM_MIXING_PDM) {
                        nummix = max;
                        numslide = nummix + 1;
                    }
                    else{
                        nummix = GR_Pulay::NUM_MIXING_PDM;
                        numslide = nummix;
                    }

                    // alpha from RDM
                    Eigen::MatrixXd a(nummix, nummix);

                    for (auto scfi = 0; scfi <= nummix; scfi++) {
                        for (auto scfj = scfi; scfj <= nummix; scfj++) {
                            a(scfi, scfj) = 0.0;

                            auto const size = ryarray_[scfi].size();
                            for (auto i = 0; i < size; i++){
                                auto const sum1 = ryarray_[scfi][i];
                                auto const sum2 = ryarray_[scfj][i];
                                a(scfi, scfj) += sum1 * sum2;
                            }
                            
                            a(scfj, scfi) = a(scfi, scfj);
                        }
                    }

					Eigen::MatrixXd ia(a.inverse());

                    auto denominator = 0.0;
                    for (auto scfi = 0; scfi <= nummix; scfi++) {
                        for (auto scfj = scfi; scfj <= nummix; scfj++) {
                            denominator += ia(scfi, scfj);
                        }
                    }

                    auto sum = 0.0;
                    std::vector<double> alden(nummix + 1);
                    for (auto scfi = 0; scfi <= nummix; scfi++){
                        auto numerator = 0.0;
                        for (auto scfj = 0; scfj <= nummix; scfj++){
                            numerator += ia(scfj, scfi);
                        }
                        alden[scfi] = numerator / denominator;
                        sum += alden[scfi];
                    }

                    // Calculate an optimized residual rho
                    auto const xmesh_size = x.size();
                    std::vector<double> optry(size);

                    auto optnorm_ry = 0.0;

                    for (auto i = 0U; i < xmesh_size; i++) {
                        optry[i] = 0.0;
                        
                        double dx;
                        if (!i) {
                            dx = x[0];
                        }

                        dx = x[i + 1] - x[i];

                        for (auto pscfiter = 0; pscfiter <= nummix; pscfiter++) {
                            if (!pscfiter) {
                                optry[i] += alden[pscfiter] * ryarray_[pscfiter][i];
                            }

                            optry[i] += alden[pscfiter] * ryarray_[pscfiter][i];
                        }

                        optnorm_ry += optry[i] * optry[i] * dx;
                    }

                    double coef_optrdm = 0.0;
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
                        for (auto pscf_iter = 0; pscf_iter <= nummix; pscf_iter++){
                            if (!pscf_iter) {
                                yarray_[0][i] += alden[pscf_iter] * py[i];
                            }

                            yarray_[0][i] += alden[pscf_iter] * yarray_[pscf_iter][i];
                        }

                        // Correction by the optimized rho
                        yarray_[0][i] += coef_optrdm * optry[i];
                        if (yarray_[0][i] < 0.0) {
                            yarray_[0][i] = 1.0e-15;
                        }
                    }

                    // Shift of rho
                    for (auto pscf_iter = numslide; 0 < pscf_iter; pscf_iter--){
                        yarray_[pscf_iter] = yarray_[pscf_iter - 1];
                    }

                    // Shift of residual rho
                    for (auto pscf_iter = numslide; 1 < pscf_iter; pscf_iter--){
                        ryarray_[pscf_iter] = ryarray_[pscf_iter - 1];
                    }
                }

				return yarray_[0];

                break;
            }
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

