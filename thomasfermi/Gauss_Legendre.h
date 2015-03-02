#include "mkl_allocator.h"
#include <vector>
#include <functional>
#include <cstdint>
#include <intrin.h>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace thomasfermi {
	namespace FEM_ALL {
		bool availableAVX();

		class Gauss_Legendre
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			Gauss_Legendre(const Gauss_Legendre &) = delete;
			Gauss_Legendre & operator=(const Gauss_Legendre &) = delete;
			Gauss_Legendre() = delete;
#endif

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
			static constexpr double EPS = 1.0E-15;
#else
			static const double EPS;
#endif
			typedef std::vector<double, mkl_allocator<double>> dvector;

			const std::size_t n_;
			const bool avxSupported_;
			dvector x_;
			dvector w_;
			void gauleg(bool usecilk);

		public:
			explicit Gauss_Legendre(std::size_t n, bool usecilk = false);
			double qgauss(double x1, double x2, bool useSSEorAVX, const std::function<double (double)> & func) const;
			const dvector & getx() const
			{ return x_; }
			const dvector & getw() const
			{ return w_; }
		};

		inline bool availableAVX()
		{
#if (_MSC_FULL_VER >= 160040219)
			std::int32_t cpuInfo[4];
			__cpuid(cpuInfo, 1);
 
			const bool osUsesXSAVE_XRSTORE = cpuInfo[2] & (1 << 27) || false;
			const bool cpuAVXSuport = cpuInfo[2] & (1 << 28) || false;
 
			if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
			{
				// Check if the OS will save the YMM registers
				const std::uint64_t xcrFeatureMask = _xgetbv(_XCR_XFEATURE_ENABLED_MASK);
				return (xcrFeatureMask & 0x6) || false;
			}
#endif
			return false;
		}

		inline Gauss_Legendre::Gauss_Legendre(std::size_t n, bool usecilk)
			:	n_(n), avxSupported_(availableAVX()),
				x_(n), w_(n)
		{
			gauleg(usecilk);
		}
	}
}
