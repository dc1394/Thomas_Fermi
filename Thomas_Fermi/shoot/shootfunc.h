#ifndef _SHOOTFUNC_H_
#define _SHOOTFUNC_H_

#ifdef _MSC_VER
	#pragma once
#endif

#pragma warning(disable : 4819)
#define _SCL_SECURE_NO_WARNINGS

#include <array>
#include <boost/numeric/ublas/matrix.hpp>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace thomasfermi {
	namespace shoot {
		class shootfunc
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			shootfunc(const shootfunc &) = delete;
			shootfunc & operator=(const shootfunc &) = delete;
			shootfunc() = delete;
#endif
		public:
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
			static constexpr std::size_t N2 = 1;
			static constexpr std::size_t NVAR = 2;
#else
			static const std::size_t N2 = 1;
			static const std::size_t NVAR = 2;
#endif
			typedef boost::numeric::ublas::vector<double> dblasvector;
			typedef std::array<double, N2> tmpary;
			typedef std::array<double, NVAR> state_type;

			static const shootfunc::tmpary V1;
			static const shootfunc::tmpary DELV;

			static shootfunc::state_type load1(double x1, const shootfunc::tmpary & v1);
			static shootfunc::dblasvector score(const shootfunc::state_type & y);
			static void rhs(const shootfunc::state_type & y, shootfunc::state_type & dydx, const double x);
		};
	}
}

#endif	// _SHOOTFUNC_H_
