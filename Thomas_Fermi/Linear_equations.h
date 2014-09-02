#pragma warning(disable : 4819)

#include "mkl_allocator.h"
#include <tuple>
#include <vector>

#if !defined(__INTEL_COMPILER) || !defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace Thomas_Fermi {
	namespace FEM_ALL {
		class Linear_equations
#if !defined(__INTEL_COMPILER) || !defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			Linear_equations & operator=(const Linear_equations &) = delete;
			Linear_equations() = delete;
#endif
			typedef std::vector<std::size_t> sivector;
			typedef std::vector<double, mkl_allocator<double>> dvector;

			const std::size_t n_;
			dvector a1_;
			const dvector a1back_;
			dvector a2_;
			const dvector a2back_;
			dvector b_;

		public:
			explicit Linear_equations(const std::tuple<dvector, dvector, dvector> & res)
				:	n_(std::get<0>(res).size()), a1_(std::get<0>(res)),
					a1back_(a1_), a2_(std::get<1>(res)),
					a2back_(a2_), b_(std::get<2>(res)) {}
			void reset(const dvector & b);
			void bound(std::size_t n_bc_given,						// 既知量の数
					   const sivector & i_bc_given,					// 既知量のインデックス
					   std::size_t n_bc_nonzero,					// 非零の既知量の数
					   const sivector & i_bc_nonzero,				// 非零の既知量のインデックス
					   const std::vector<double> & v_bc_nonzero);	// 非零の既知量
			Linear_equations::dvector LEsolver();
		};
	}
}
