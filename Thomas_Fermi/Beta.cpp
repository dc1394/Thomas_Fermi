#include "Beta.h"
#include <stdexcept>
#include <cstdint>

namespace Thomas_Fermi {
	namespace FEM_ALL {
		double Beta::operator()(double x) const
		{
			std::uint32_t klo = 0;												// 表の中の正しい位置を二分探索で求める
			std::uint32_t khi = static_cast<std::uint32_t>(size_ - 1);

			while (khi - klo > 1) {
				const std::uint32_t k = static_cast<std::uint32_t>((khi + klo) >> 1);

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
