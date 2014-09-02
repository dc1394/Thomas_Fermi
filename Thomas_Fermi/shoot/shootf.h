#pragma warning(disable : 4819)
// Visual C++ 2008以下は対応しません、GCCはC++11オプションで動かして下さい
#if (_MSC_VER <= 1500) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
	#include <boost/static_assert.hpp>
	BOOST_STATIC_ASSERT(false);
#endif

#include "load2.h"
#include <tuple>
#include <utility>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace Thomas_Fermi {
	namespace shoot {
		class shootf
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			shootf(const shootf &) = delete;
			shootf & operator=(const shootf &) = delete;
			shootf() = delete;
#endif

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
			static constexpr double EPS = 1.0E-14;
#else
			static const double EPS;
#endif
			typedef std::vector<double> dvector;

		public:
			typedef std::tuple<dvector, const dvector, const shootfunc::tmpary> result_type;

		private:
			typedef std::function<shootfunc::state_type (double, const shootfunc::tmpary &)> loadfunctype;
			typedef std::function<shootfunc::dblasvector (const shootfunc::state_type & )> scorefunctype;

			shootfunc::tmpary v1_;
			shootfunc::tmpary v2_;
			const shootfunc::tmpary delv1_;
			const shootfunc::tmpary delv2_;

			const double dx_;
			const double eps_;

			const loadfunctype load1_;
			const loadfunctype load2_;
			const scorefunctype score_;

			shootf::result_type createResult(double x1, double xf,
											 const dvector & res1, const dvector & res2) const;

		public:
			shootf(const shootfunc::tmpary & v1, const shootfunc::tmpary & v2,
				   const shootfunc::tmpary & delv1, const shootfunc::tmpary & delv2,
				   double dx, double eps,
				   const loadfunctype & load1,
				   const load2 & l2,
				   const scorefunctype & score);
			shootf::result_type operator()(double x1, double x2, double xf);
		};

		inline shootf::shootf(const shootfunc::tmpary & v1, const shootfunc::tmpary & v2,
							  const shootfunc::tmpary & delv1, const shootfunc::tmpary & delv2,
							  double dx, double eps,
							  const loadfunctype & load1,
							  const load2 & l2,
							  const scorefunctype & score)
		 :	v1_(v1), v2_(v2), delv1_(delv1), delv2_(delv2),
			dx_(dx), eps_(eps), load1_(load1),
			load2_(std::bind(&load2::operator(), std::ref(l2), std::placeholders::_1, std::placeholders::_2)),
			score_(score) {}
	}
}
