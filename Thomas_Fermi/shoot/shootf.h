/*! \file shootf.h
    \brief 狙い撃ち法により、y(x)を求める

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/

#ifndef _SHOOTF_H_
#define _SHOOTF_H

#pragma once

#include "load2.h"
#include "../myfunctional/Functional.h"
#include <functional>
#include <tuple>
#include <utility>
#include <boost/numeric/ublas/matrix.hpp>

namespace thomasfermi {
	namespace shoot {
		class shootf final
		{
			shootf(const shootf &) = delete;
			shootf & operator=(const shootf &) = delete;
			shootf() = delete;

			static constexpr double EPS = 1.0E-14;

			typedef std::vector<double> dvector;

		public:
			typedef std::tuple<dvector, const dvector, const shootfunc::tmpary> result_type;

		private:
            typedef myfunctional::Functional<shootfunc::state_type(double, shootfunc::tmpary const &)> loadfunctype;
            typedef myfunctional::Functional<shootfunc::dblasvector(const shootfunc::state_type &)> scorefunctype;

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

#endif  // _SHOOTF_H_
