#include "FOElement.h"

namespace thomasfermi {
	namespace FEM_ALL {
		FOElement::FOElement(std::size_t nint, bool useSSEorAVX, bool usecilk, const dvector & coords, dvector && beta)
			:	FEM(nint, useSSEorAVX, usecilk, coords, std::move(beta)),
				N1_([](double r) { return 0.5 * (1.0 - r); }),
				N2_([](double r) { return 0.5 * (1.0 + r); }),
				fun1_([this](double r, double xl, std::size_t ielem)
				{ return - N1_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]]) * xl * 0.5; }),
				fun2_([this](double r, double xl, std::size_t ielem)
				{ return - N2_(r) * func_(N1_(r) * coords_[lnods_[0][ielem]] + N2_(r) * coords_[lnods_[1][ielem]]) * xl * 0.5; })
		{
			ntnoel_ = 2;
			nelem_ = nnode_ - 1;
			
			initialize();
			
			for (std::size_t i = 0; i < nelem_; i++) {
				lnods_[0][i] = i;
				lnods_[1][i] = i + 1;
			}
		}

		FEM::dvector FOElement::getdndr() const
		{
			dvector dndr(ntnoel_);
			dndr[0] = - 0.5;
			dndr[1] = 0.5;

			return std::move(dndr);
		}

		FEM::dvector FOElement::getc(std::size_t ielem) const
		{
			const double xl = coords_[lnods_[1][ielem]] - coords_[lnods_[0][ielem]];
			dvector c(ntnoel_);

			c[0] = gl_.qgauss(- 1.0, 1.0, useSSEorAVX_, [this, xl, ielem](double r){ return fun1_(r, xl, ielem); });
			c[1] = gl_.qgauss(- 1.0, 1.0, useSSEorAVX_, [this, xl, ielem](double r){ return fun2_(r, xl, ielem); });

			return std::move(c);
		}
	}
}
