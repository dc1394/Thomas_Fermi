/*! \file foelement.cpp
	\brief 一次要素のクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	(but this is originally adapted by 渡辺浩志 for stiff5.c from http://www.sml.k.u-tokyo.ac.jp/members/nabe/FEM/FEM.pdf )
	This software is released under the BSD-2 License.
*/

#include "foelement.h"
#include "myfunctional/functional.h"

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ

		FOElement::FOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usesimd, bool usetbb) :
			FEM(std::move(beta), coords, nint, usesimd, usetbb),
			fun1_([this](double r, double xl, std::size_t ielem)
			{ return -N1_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]]) * xl * 0.5; }),
			fun2_([this](double r, double xl, std::size_t ielem)
			{ return -N2_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]]) * xl * 0.5; }),
			N1_([](double r) { return 0.5 * (1.0 - r); }),
			N2_([](double r) { return 0.5 * (1.0 + r); })
		{
			ntnoel_ = 2;
			nelem_ = nnode_ - 1;
			
			initialize();
			
			for (auto i = 0U; i < nelem_; i++) {
				lnods_[0][i] = i;
				lnods_[1][i] = i + 1;
			}
		}

		// #endregion コンストラクタ

		// #region メンバ関数

		FEM::dvector FOElement::getc(std::size_t ielem) const
		{
			const double xl = coords_[lnods_[1][ielem]] - coords_[lnods_[0][ielem]];
			dvector c(ntnoel_);

			c[0] = gl_.qgauss(
				myfunctional::make_functional([this, xl, ielem](double r){ return fun1_(r, xl, ielem); }),
				usesimd_,
				-1.0,
				1.0);

			c[1] = gl_.qgauss(
				myfunctional::make_functional([this, xl, ielem](double r){ return fun2_(r, xl, ielem); }),
				usesimd_,
				-1.0,
				1.0);

			return std::move(c);
		}

		FEM::dvector FOElement::getdndr() const
		{
			dvector dndr(ntnoel_);
			dndr[0] = - 0.5;
			dndr[1] = 0.5;

			return std::move(dndr);
		}

		// #endregion メンバ関数
	}
}
