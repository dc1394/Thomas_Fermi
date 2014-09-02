#include "../Beta.h"
#include <tuple>
#include <memory>
#include <fstream>
#include <utility>
#include <cstdint>
#include <boost/cast.hpp>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace Thomas_Fermi {
	namespace MakeRhoEn {
		class MakeRhoEnergy 
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
			: private boost::noncopyable
#endif
		{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
			MakeRhoEnergy(const MakeRhoEnergy &) = delete;
			MakeRhoEnergy & operator=(const MakeRhoEnergy &) = delete;
			MakeRhoEnergy() = delete;
#endif
			typedef std::tuple<std::vector<double>, std::shared_ptr<const FEM_ALL::Beta>,
							   double> parameter_type;
			static const double ALPHA;

			const std::vector<double> xvec_;
			const std::shared_ptr<const FEM_ALL::Beta> pbeta_;
			
			std::ofstream ofs;
			
			const std::size_t size_;

			const double v1_;
			const double dx_;

			const std::int32_t max_;

			double a;

			double y(double x) const;
			double rhofunc(double x) const;
			double makeEnergy() const;
			void saveresult1();
			void saveresult2(std::size_t n, bool useAVXorSSE, bool usecilk);
			void saveresult3();

		public:
			MakeRhoEnergy::MakeRhoEnergy(const parameter_type & pt);
			void saveresult(std::size_t n, bool useAVXorSSE, bool usecilk);
		};

		inline MakeRhoEnergy::MakeRhoEnergy(const parameter_type & pt)
			:	xvec_(std::get<0>(pt)), pbeta_(std::get<1>(pt)),
				size_(xvec_.size()), v1_(std::get<2>(pt)),
				dx_((xvec_[2] - xvec_[1]) * 2.0),
				max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] /
					 MakeRhoEnergy::ALPHA / dx_))
		{
		}

#ifdef __GXX_EXPERIMENTAL_CXX0X__
		constexpr
#endif
		double power(double x, std::uint32_t y);
	}
}
