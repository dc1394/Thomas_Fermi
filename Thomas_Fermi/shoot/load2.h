#include "shootfunc.h"
#include "../Spline.h"
#include <boost/optional.hpp>

namespace Thomas_Fermi {
	namespace shoot {
		class load2	{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
			static constexpr std::size_t XAYSIZE = 150;
			static constexpr double THRESHOLD = 65.0;
			static constexpr double k = 0.190785707092222;
			static constexpr double alpha = 3.886;
#else
			static const std::size_t XAYSIZE = 150;
			static const double THRESHOLD;
			static const double k;
			static const double alpha;
#endif
			boost::optional<Spline> psplint;

			static double phiT(double x);
			static double dphiT(double x);

		public:
			load2();
			~load2() { psplint = boost::none; }
			shootfunc::tmpary make_v2(double x2) const;
			shootfunc::state_type operator()(double x2, const shootfunc::tmpary & v2) const;
		};
	}
}
