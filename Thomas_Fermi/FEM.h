﻿#ifndef _FEM_H_
#define _FEM_H_

#pragma once

#ifdef _MSC_VER
    #pragma warning(disable : 4819)
    #define _SCL_SECURE_NO_WARNINGS
#endif

#include "Beta.h"
#include "gausslegendre/Gauss_Legendre.h"
#include "mkl_allocator.h"
#include <functional>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/optional.hpp>

namespace thomasfermi {
	namespace FEM_ALL {
		class FEM final
		{
		public:
			typedef std::vector<double> dvector;
			typedef std::vector<double, mkl_allocator<double>> dmklvector;

		protected:
			const std::size_t nnode_;
			const bool useSSEorAVX_;
			const bool usecilk_;

		private:
			typedef boost::multi_array<double, 2> dmatrix;
			const std::size_t nint_;

		protected:
			std::size_t ntnoel_;
			std::size_t nelem_;
			boost::multi_array<std::size_t, 2> lnods_;				

		private:
			dmklvector a1_;
			dmklvector a2_;
			dmatrix astiff_;
			dmklvector b_;

		protected:
			const dvector coords_;

		private:
			const dvector beta_;

		protected:
			Gauss_Legendre gl_;
			std::shared_ptr<const Beta> pbeta_;
			std::function<double (double)> func_;
			void initialize();

		private:
			virtual dvector getdndr() const = 0;
			virtual dvector getc(std::size_t ielem) const = 0;

			void element(std::size_t ielem);
			void amerge(std::size_t ielem);

		public:
			FEM(std::size_t nint, bool useSSEorAVX, bool usecilk, const dvector & coords, dvector && beta);
			virtual ~FEM() {}
			void stiff();
			void reset(const dvector & beta);
			void stiff2();

			std::size_t getntnoel() const
			{ return ntnoel_; }

			std::size_t getnnode() const
			{ return nnode_; }

			std::tuple<dmklvector, dmklvector, dmklvector> createresult() const
			{ return std::make_tuple(a1_, a2_, b_); };

			const std::shared_ptr<const Beta> & getpbeta() const
			{ return pbeta_; }

			const FEM::dmklvector & getb() const
			{ return b_; }

            FEM(const FEM &) = delete;
            FEM & operator=(const FEM &) = delete;
            FEM() = delete;
		};

		inline FEM::FEM(std::size_t nint, bool useSSEorAVX, bool usecilk, const dvector & coords, dvector && beta)
			:	nnode_(coords.size()),
				useSSEorAVX_(useSSEorAVX),
				usecilk_(usecilk),
				nint_(nint),
				a1_(nnode_, 0.0),
				a2_(nnode_ - 1, 0.0),
				b_(nnode_, 0.0),
				coords_(coords),
				beta_(std::move(beta)),
				gl_(nint_, usecilk_),
				pbeta_(std::make_shared<const Beta>(coords_, beta_)),
				func_(std::bind(&Beta::operator(), std::ref(*pbeta_), std::placeholders::_1))
		{
			BOOST_ASSERT(coords_.size() == beta_.size());
		}
	}
}

#endif  // _FEM_H_