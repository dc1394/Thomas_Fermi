#include "SOElement.h"

namespace Thomas_Fermi {
	namespace FEM_ALL {
		SOElement::SOElement(std::size_t nint, const dvector & coords, dvector && beta)
			:	FEM(nint, coords, std::move(beta)),
				N1_([](double r) { return - 0.5 * r * (1.0 - r); }),
				N2_([](double r) { return 0.5 * r * (1.0 + r); }),
				N3_([](double r) { return (1.0 - r * r); })
		{
			ntnoel_ = 3;
			nelem_ = (nnode_ - 1) >> 1;
			
			initialize();

			for (std::size_t i = 0; i < nelem_; i++) {
				lnods_[0][i] = 2 * i;
				lnods_[1][i] = 2 * i + 2;
				lnods_[2][i] = 2 * i + 1;
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
			const double xl = coords_[lnods_[1][ielem]] - coords_[lnods_[0][ielem]];

			dvector c(ntnoel_);
			c[0] = gl_.qgauss(- 1.0, 1.0, [this, ielem](double r)
				   { return - N1_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]] + N3_(r) * coords_[lnods_[2][ielem]]); },
				   true) * xl * 0.5;
			c[1] = gl_.qgauss(- 1.0, 1.0, [this, ielem](double r)
				   { return - N2_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]] + N3_(r) * coords_[lnods_[2][ielem]]); },
				   true) * xl * 0.5;
			c[0] = gl_.qgauss(- 1.0, 1.0,  [this, ielem](double r)
				   { return - N3_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]] + N3_(r) * coords_[lnods_[2][ielem]]); },
				   true) * xl * 0.5;

			return c;
		}
	}
}
