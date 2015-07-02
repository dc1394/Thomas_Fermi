/*! \file soelement.h
	\brief 二次要素のクラスの実装

	Copyright ©  2015 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/


#include "SOElement.h"

namespace thomasfermi {
	namespace femall {
		// #region コンストラクタ

		SOElement::SOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk)
			:	FEM(std::move(beta), coords, nint, usecilk),
				a2_(nnode_ - 2, 0.0)
		{
			auto const N1tmp = [](double r) { return -0.5 * r * (1.0 - r); };
			N1_ = std::cref(N1tmp);

			auto const N2tmp = [](double r) { return 0.5 * r * (1.0 + r); };
			N2_ = std::cref(N2tmp);
			
			auto const N3tmp = [](double r) { return (1.0 - r * r); };
			N3_ = std::cref(N3tmp);

			ntnoel_ = 3;
			nelem_ = (nnode_ - 1) >> 1;
			
			initialize();

			for (auto i = 0U; i < nelem_; i++) {
                (*plnods_)[i][0] = 2 * i;
                (*plnods_)[i][1] = 2 * i + 2;
                (*plnods_)[i][2] = 2 * i + 1;
			}
		}
		
		// #endregion コンストラクタ

		// #region publicメンバ関数

		std::tuple<FEM::dmklvector, FEM::dmklvector, FEM::dmklvector, FEM::dmklvector> SOElement::createresult() const
		{
			//FEM::resultmap mymap;

			//mymap["a0"] = a0_;
			//mymap["a1"] = a1_;
			//mymap["a2"] = a2_;
			//mymap["b"] = b_;

			return std::make_tuple(a0_, a1_, a2_, b_);
		}
		
		// #endregion publicメンバ関数

		// #region privateメンバ関数

		void SOElement::amerge(std::size_t ielem)
		{
            for (auto i = 0UL; i < ntnoel_; i++) {
                for (auto j = 0UL; j < ntnoel_; j++) {
                    auto const lnodi = (*plnods_)[ielem][i];
                    auto const lnodj = (*plnods_)[ielem][j];
                    if (lnodi == lnodj) {
                        a0_[lnodi] += astiff_[i][j];
                    }
                    else if (lnodi == lnodj - 1 && lnodi < nnode_ - 1) {
                        a1_[lnodi] += astiff_[i][j];
                    }
                    else if (lnodi == lnodj - 2 && lnodi < nnode_ - 2) {
                        a2_[lnodi] += astiff_[i][j];
                    }
                }
            }
		}
		
		void SOElement::element(std::size_t ielem)
		{
			astiffclear();

			for (auto ir = 0U; ir < nint_; ir++) {
				auto const dndr(getdndr(gl_.X()[ir]));

				FEM::element(dndr, ielem, ir);
			}
		}

		FEM::dvector SOElement::getdndr(double r) const
		{
			dvector dndr(ntnoel_);
			dndr[0] = r - 0.5;
			dndr[1] = r + 0.5;
			dndr[2] = - 2.0 * r;

			return dndr;
		}

		FEM::dvector SOElement::getc(std::size_t ielem) const
		{
            auto const xl = coords_[(*plnods_)[ielem][1]] - coords_[(*plnods_)[ielem][0]];

			dvector c(ntnoel_);
			c[0] = gl_.qgauss(
				myfunctional::make_functional(
					[this, ielem](double r)
				   { return - N1_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]] + N3_(r) * coords_[(*plnods_)[ielem][2]]); }),
				   -1.0,
				   1.0) * xl * 0.5;
			c[1] = gl_.qgauss(
				myfunctional::make_functional([this, ielem](double r)
				   { return - N2_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]] + N3_(r) * coords_[(*plnods_)[ielem][2]]); }),
				   -1.0,
				   1.0) * xl * 0.5;
			c[2] = gl_.qgauss(
				myfunctional::make_functional([this, ielem](double r)
				   { return - N3_(r) * func_(N1_(r) * coords_[(*plnods_)[ielem][0]] + N2_(r) * coords_[(*plnods_)[ielem][1]] + N3_(r) * coords_[(*plnods_)[ielem][2]]); }),
				   -1.0,
				   1.0) * xl * 0.5;

			return c;
		}

		// #endregion privateメンバ関数
	}
}
