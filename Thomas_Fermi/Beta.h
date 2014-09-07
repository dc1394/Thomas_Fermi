#ifndef _BETA_H_
#define _BETA_H_

#ifdef _MSC_VER
	#pragma once
#endif

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

#include "mkl_allocator.h"
#include <vector>

namespace thomasfermi {
	namespace FEM_ALL {
		class Beta
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			Beta(const Beta &) = delete;
			Beta & operator=(const Beta &) = delete;
			Beta() = delete;
#endif
			const std::size_t size_;
			const std::vector<double> xvec_;
			const std::vector<double> yvec_;

		public:
			Beta(const std::vector<double> & xvec, const std::vector<double> & yvec)
				:	size_(xvec.size()), xvec_(xvec), yvec_(yvec) {}
			double operator()(double x) const;
		};
	}
}

#endif	// _BETA_H_
