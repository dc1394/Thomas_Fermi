#include "FEM.h"
#include <cstdint>
#include <boost/cast.hpp>
#include <boost/assert.hpp>

#if !defined(__INTEL_COMPILER) && defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/foreach.hpp>
#endif

#ifdef __cilk
	#include <cilk/cilk.h>
#else
	#define cilk_for for
#endif

namespace Thomas_Fermi {
	namespace FEM_ALL {
		void FEM::initialize()
		{
			astiff_.resize(boost::extents[ntnoel_][ntnoel_]);
			lnods_.resize(boost::extents[ntnoel_][nelem_]);
		}

		void FEM::element(std::size_t ielem)
		{
			for (std::size_t i = 0; i < ntnoel_; i++)
				for (std::size_t j = 0; j < ntnoel_; j++)
					astiff_[i][j] = 0.0;

			for (std::size_t ir = 0; ir < nint_; ir++) {
				const dvector dndr(getdndr());
				double ajacob = 0.0;

				for (std::size_t i = 0; i < ntnoel_; i++)
					ajacob += dndr[i] * coords_[lnods_[i][ielem]];

				const double detjac = ajacob;
				const double ajainv = 1.0 / ajacob;
				
				dvector dndx(ntnoel_);
				for (std::size_t i = 0; i < ntnoel_; i++)
					dndx[i] = dndr[i] * ajainv;

				const double detwei = detjac * gl_.getw()[ir];
				for (std::size_t i = 0; i < ntnoel_; i++)
					for (std::size_t j = 0; j < ntnoel_; j++)
						astiff_[i][j] += detwei * dndx[i] * dndx[j];
			}
		}

		void FEM::amerge(std::size_t ielem)
		{
			a1_[ielem] += astiff_[0][0];
			a1_[ielem + 1] += astiff_[1][1];
			a2_[ielem] = astiff_[0][1];
		}

		void FEM::stiff()
		{
			const std::int32_t size = boost::numeric_cast<std::int32_t>(nelem_);
			for (std::int32_t ielem = 0; ielem < size; ielem++)
				element(ielem);

			if (usecilk_) {
				cilk_for (std::int32_t ielem = 0; ielem < size; ielem++) {
					amerge(ielem);

					const dvector c(getc(ielem));
					for (std::size_t i = 0; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				}
			} else {
				for (std::int32_t ielem = 0; ielem < size; ielem++) {
					amerge(ielem);

					const dvector c(getc(ielem));
					for (std::size_t i = 0; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				}
			}
		}

		void FEM::reset(const dvector & beta)
		{
			pbeta_.reset();
			pbeta_ = std::make_shared<const Beta>(coords_, beta);
			func_ = std::bind(&Beta::operator(), std::ref(*pbeta_), std::placeholders::_1);

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
#ifdef __INTEL_COMPILER
#pragma loop count min(1024)
#endif
			for (double & v : b_)
				v = 0.0;
#else
			BOOST_FOREACH (double & v, b_)
				v = 0.0
#endif
		}

		void FEM::stiff2()
		{
			const std::int32_t size = boost::numeric_cast<std::int32_t>(nelem_);
			if (usecilk_) {
				cilk_for (std::int32_t ielem = 0; ielem < size; ielem++) {
					const dvector c(getc(ielem));
					for (std::size_t i = 0; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				}
			} else {
				for (std::int32_t ielem = 0; ielem < size; ielem++) {
					const dvector c(getc(ielem));
					for (std::size_t i = 0; i < ntnoel_; i++)
						b_[lnods_[i][ielem]] += c[i];
				}
			}
		}
	}
}
