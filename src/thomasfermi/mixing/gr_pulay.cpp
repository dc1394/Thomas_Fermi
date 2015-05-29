/*! \file gr_pulay.cpp
	\brief guaranteed-reduction Pulay法でyの合成を行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
    (but this is originally adapted by T.Ozaki for GR_Pulay.c from ADPACK)
	This software is released under the BSD-2 License.
*/

#include "gr_pulay.h"

namespace thomasfermi {
    namespace mixing {
        // #region コンストラクタ

        GR_Pulay::GR_Pulay(std::shared_ptr<Data> const & pdata) :
            SimpleMixing(pdata)
        {
        }

        // #region コンストラクタ

        // #region publicメンバ関数

        femall::FEM::dmklvector GR_Pulay::operator()(std::int32_t scfiter, femall::FEM::dmklvector const & y)
        {
            auto const size = y.size();
            BOOST_ASSERT(size == Ybefore.size());

            switch (scfiter) {
            case 0:
                setyryarray(y);
                return yarray[0];
                break;

            case 1:
                ryarray[2] = getry();

                yarray[0] = SimpleMixing::operator()(yarray[0], yarray[1]);
                yarray[2] = yarray[1];
                yarray[1] = yarray[0];
                
                return yarray[0];

            default:
                break;
            }
            femall::FEM::dmklvector newy(size);

            for (auto i = 0U; i < size; i++) {
                newy[i] = Ybefore()[i] + pdata_->iteration_mixing_weight_ * (y[i] - Ybefore()[i]);
            }
            
            Ybefore = y;

            return newy;
        }

        // #endregion publicメンバ関数
        
        // #region privateメンバ関数

        std::vector<double> GR_Pulay::getry()
        {
            auto const size = yarray[0].size();
            BOOST_ASSERT(size == yarray[1].size());

            std::vector<double> ry(size);
            for (auto i = 0U; i < size; i++) {
                ry[i] = yarray[0][i] - yarray[1][i];
            }

            return ry;
        }

        std::vector<double> GR_Pulay::getry(femall::FEM::dmklvector const & y)
        {
            auto const size = y.size();
            BOOST_ASSERT(size == Ybefore.size());

            std::vector<double> ry(size);
            for (auto i = 0U; i < size; i++) {
                ry[i] = y[i] - Ybefore()[i];
            }

            return ry;
        }

        void GR_Pulay::setyryarray(femall::FEM::dmklvector const & y)
        {
            ryarray[2] = getry(y);
            yarray[0] = SimpleMixing::operator()(y);
            yarray[1] = y;
            yarray[2] = Ybefore;
        }

        // #endregion privateメンバ関数
    }
}

