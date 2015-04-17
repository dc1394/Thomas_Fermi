
#ifdef _MSC_VER
#pragma warning(disable : 4819)
#define _SCL_SECURE_NO_WARNINGS
#endif


#include "mkl_allocator.h"
#include <tuple>
#include <vector>

namespace thomasfermi {
	namespace fem_all {
		class Linear_equations final {
			Linear_equations & operator=(const Linear_equations &) = delete;
			Linear_equations() = delete;
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
