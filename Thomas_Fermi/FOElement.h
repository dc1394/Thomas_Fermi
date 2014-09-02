#include "FEM.h"
#include <functional>

namespace Thomas_Fermi {
	namespace FEM_ALL {
		class FOElement : public FEM {
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			FOElement(const FOElement &) = delete;
			FOElement & operator=(const FOElement &) = delete;
			FOElement() = delete;
#endif
			const std::function<double (double)> N1_;
			const std::function<double (double)> N2_;
			const std::function<double (double, double, std::size_t)> fun1_;
			const std::function<double (double, double, std::size_t)> fun2_;

			dvector getdndr() const override;
			dvector getc(std::size_t ielem) const override;

		public:
			FOElement(std::size_t nint, bool useSSEorAVX, bool usecilk, const dvector & coords, dvector && beta);
			virtual ~FOElement() {}
		};
	}
}
