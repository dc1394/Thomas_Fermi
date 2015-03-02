/*! \file Beta.h
\brief β(x)を計算するクラスの実装

Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#include "Beta.h"
#include <cstdint>
#include <stdexcept>

namespace thomasfermi {
	namespace FEM_ALL {
		double Beta::operator()(double x) const
		{
			auto klo = 0;
			auto khi = static_cast<std::uint32_t>(size_ - 1);

            // 表の中の正しい位置を二分探索で求める
			while (khi - klo > 1) {
				auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

				if (xvec_[k] > x)
					khi = k;
				else
					klo = k;
			}

			// yvec[i] = f(xvec[i]), yvec[i + 1] = f(xvec[i + 1])の二点を通る直線を代入
			return (yvec_[khi] - yvec_[klo]) / (xvec_[khi] - xvec_[klo]) * (x - xvec_[klo]) + yvec_[klo];
		}
	}
}
