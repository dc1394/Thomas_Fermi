/*! \file simplemixing.cpp
	\brief 一次混合法でyの合成を行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "simplemixing.h"

namespace thomasfermi {
    namespace mixing {
        SimpleMixing::SimpleMixing(std::shared_ptr<Data> const & pdata) :
            Yold(
				[this] { return std::cref(yold_); },
                [this](femall::FEM::dmklvector const & val) { 
                    yold_ = val;
                    return val;
            }),
            pdata_(pdata)
        {
        }

		femall::FEM::dmklvector SimpleMixing::operator()(femall::FEM::dmklvector const & y)
        {
            auto const size = y.size();
            BOOST_ASSERT(size == yold_.size());

            femall::FEM::dmklvector newy(size);

			for (auto i = 0U; i < size; i++) {
                newy[i] = yold_[i] + pdata_->iteration_mixing_weight_ * (y[i] - yold_[i]);
            }
            
            yold_ = y;

            return newy;
        }
    }
}

