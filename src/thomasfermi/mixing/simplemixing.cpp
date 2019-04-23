/*! \file simplemixing.cpp
	\brief 一次混合法でyの合成を行うクラスの実装
	Copyright © 2015-2019 @dc1394 All Rights Reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "simplemixing.h"

namespace thomasfermi {
    namespace mixing {
        SimpleMixing::SimpleMixing(std::shared_ptr<Data> const & pdata)
            :   Yold(
				[this] { return std::cref(yold_); },
                [this](std::vector<double> const & val) { 
                    yold_ = val;
                    return val;
                }),
                pdata_(pdata)
        {
        }

		std::vector<double> SimpleMixing::operator()(std::vector<double> const & y)
        {
            auto const size = y.size();
            BOOST_ASSERT(size == yold_.size());

            std::vector<double> newy(size);

			for (auto i = 0U; i < size; i++) {
                newy[i] = yold_[i] + pdata_->iteration_mixing_weight_ * (y[i] - yold_[i]);
            }
            
            yold_ = y;

            return newy;
        }
    }
}
