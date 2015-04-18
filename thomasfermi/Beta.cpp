/*! \file beta.cpp
	\brief β(x)を計算するクラスの実装

	Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#include "beta.h"
#include <cstdint>	// for std::uint32_t

namespace thomasfermi {
	namespace fem_all {
		// #region コンストラクタ

		Beta::Beta(std::vector<double> const & xvec, std::vector<double> const & yvec) :
			size_(xvec.size()),
			xvec_(xvec),
			yvec_(yvec)
		{
		}
		
		// #endregion コンストラクタ
		
		// #region メンバ関数

		double Beta::operator()(double x) const
		{
			auto klo = 0U;
			auto khi = static_cast<std::uint32_t>(size_ - 1);

            // 表の中の正しい位置を二分探索で求める
			while (khi - klo > 1) {
				auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

				if (xvec_[k] > x) {
					khi = k;
				}
				else {
					klo = k;
				}
			}

			// yvec[i] = f(xvec[i]), yvec[i + 1] = f(xvec[i + 1])の二点を通る直線を代入
			return (yvec_[khi] - yvec_[klo]) / (xvec_[khi] - xvec_[klo]) * (x - xvec_[klo]) + yvec_[klo];
		}

		// #endregion メンバ関数
	}
}
