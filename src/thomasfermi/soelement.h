#include "FEM.h"
#include <functional>

namespace Thomas_Fermi {
	namespace FEM_ALL {
		class SOElement : public FEM {
			std::function<double (double)> N1_;
			std::function<double (double)> N2_;
			std::function<double (double)> N3_;

			virtual dvector getdndr(double r) const;
			virtual dvector getc(std::size_t ielem) const;

		public:
			SOElement(std::size_t nint, const dvector & coords, dvector && beta);
			virtual ~SOElement() {}
		};
	}
}
