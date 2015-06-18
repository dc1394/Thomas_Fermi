#include "SOElement.h"

namespace thomasfermi {
	namespace femall {
		SOElement::SOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk)
			:	FEM(std::move(beta), coords, nint, usecilk),
				a3_(nnode_ - 1, 0.0)
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
				(*plnods_)[0][i] = 2 * i;
				(*plnods_)[1][i] = 2 * i + 2;
				(*plnods_)[2][i] = 2 * i + 1;
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
			auto const xl = coords_[(*plnods_)[1][ielem]] - coords_[(*plnods_)[0][ielem]];

			dvector c(ntnoel_);
			c[0] = gl_.qgauss(
				myfunctional::make_functional(
					[this, ielem](double r)
				   { return - N1_(r) * func_(N1_(r) * coords_[(*plnods_)[0][ielem]] + N2_(r) * coords_[(*plnods_)[1][ielem]] + N3_(r) * coords_[(*plnods_)[2][ielem]]); }),
				   -1.0,
				   1.0) * xl * 0.5;
			c[1] = gl_.qgauss(
				myfunctional::make_functional([this, ielem](double r)
				   { return - N2_(r) * func_(N1_(r) * coords_[(*plnods_)[0][ielem]] + N2_(r) * coords_[(*plnods_)[1][ielem]] + N3_(r) * coords_[(*plnods_)[2][ielem]]); }),
				   -1.0,
				   1.0) * xl * 0.5;
			c[0] = gl_.qgauss(
				myfunctional::make_functional([this, ielem](double r)
				   { return - N3_(r) * func_(N1_(r) * coords_[(*plnods_)[0][ielem]] + N2_(r) * coords_[(*plnods_)[1][ielem]] + N3_(r) * coords_[(*plnods_)[2][ielem]]); }),
				   -1.0,
				   1.0) * xl * 0.5;

			return c;
		}
	}
}
