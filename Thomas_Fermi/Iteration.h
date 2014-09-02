#include "shoot/shootfunc.h"
#include "FOElement.h"
#include "Linear_equations.h"
#include <tuple>
#include <memory>
#include <utility>
#include <boost/optional.hpp>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace Thomas_Fermi {
	namespace FEM_ALL {
		class Iteration
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			Iteration(const Iteration &) = delete;
			Iteration & operator=(const Iteration &) = delete;
			Iteration() = delete;
#endif

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
			static constexpr std::size_t N_BC_GIVEN = 2;
			static constexpr double IterationTHRESHOLD = 1.0;
			static constexpr double IterationREDUCTION = 0.15;
#else
			static const std::size_t N_BC_GIVEN = 2;
			static const double IterationTHRESHOLD;
			static const double IterationREDUCTION;
#endif
			const bool useSSEorAVX_;
			const bool usecilk_;
			const bool avxSupported_;
			const double TOL_;
			double alpha_;
			FEM::dvector x_;
			FEM::dmklvector y_;
			FEM::dmklvector ybefore_;
			std::unique_ptr<FEM> pfem_;
			std::shared_ptr<const Beta> pbeta_;
			boost::optional<Linear_equations> ple_;
			std::vector<std::size_t> i_bc_given_;
			FEM::dvector v_bc_nonzero_;

			double y1_;
			double y2_;
			double v1_;

			FEM::dvector make_beta() const;
			void ymix();
			double IterationError() const;

		public:
			typedef std::tuple<FEM::dvector, std::shared_ptr<const Beta>, double> result_type;
				// 原点に近い方、無限遠に近い方、適合点、要素の間隔
			Iteration(double x1, double x2, double xf, double dx,
				// Gauss-Legendreの積分点、ベクトル化するかどうか、並列化するかどうか
				std::size_t n, bool useSSEorAVX, bool usecilk,
				// Iterationの許容誤差、Iterationの1次混合の値α
				double TOL, double alpha);
			~Iteration() { ple_ = boost::none; }
			void Iterationloop();
			result_type makeresult()
			{ return std::make_tuple(std::move(x_), std::move(pbeta_), v1_); }
		};
	}

	template <typename T>
	inline T sqr(T x)
	{ return x * x; }
}
