/*! \file simplemixing.cpp
	\brief 一次混合法でyの合成を行うクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "simplemixing.h"

namespace thomasfermi {
    namespace mixing {
        SimpleMixing::SimpleMixing(std::shared_ptr<Data> const & pdata) :
            Ybefore(
                [this] { return ybefore_; },
                [this](femall::FEM::dmklvector const & val) { 
                    ybefore_ = val;
                    return val;
            }),
            pdata_(pdata)
        {
        }

		SimpleMixing::~SimpleMixing() noexcept
		{
		}

        femall::FEM::dmklvector SimpleMixing::operator()(femall::FEM::dmklvector const & y)
        {
            auto const size = y.size();
            BOOST_ASSERT(size == ybefore_.size());

            femall::FEM::dmklvector newy(size);

            for (auto i = 0U; i < size; i++) {
                newy[i] = ybefore_[i] + pdata_->iteration_mixing_weight_ * (y[i] - ybefore_[i]);
            }
            
            ybefore_ = y;

            return newy;
        }

        femall::FEM::dmklvector SimpleMixing::operator()(femall::FEM::dmklvector const & newy, femall::FEM::dmklvector const & oldy)
        {
            auto const size = newy.size();
            BOOST_ASSERT(size == oldy.size());

            femall::FEM::dmklvector y(size);

            for (auto i = 0U; i < size; i++) {
                y[i] = oldy[i] + pdata_->iteration_mixing_weight_ * (newy[i] - oldy[i]);
            }

            return y;
        }
    }
}

