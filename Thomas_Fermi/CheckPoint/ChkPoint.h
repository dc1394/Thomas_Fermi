#ifndef _CHKPOINT_H_
#define _CHKPOINT_H_

#include "FastArenaObject.h"
#include <memory>
#include <utility>
#include <cstdint>

#if (_MSC_VER >= 1700) || defined(__GXX_EXPERIMENTAL_CXX0X__)
    #include <chrono>
#else
    #include <boost/chrono/chrono.hpp>
#endif

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace CheckPoint {
#if (_MSC_VER >= 1700) || defined(__GXX_EXPERIMENTAL_CXX0X__)
    using namespace std::chrono;
#else
    using namespace boost::chrono;
#endif

#ifdef _MSC_VER
	void usedmem();
#endif

	class ChkPoint
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
		: private boost::noncopyable
#endif
	{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
		ChkPoint(const ChkPoint &) = delete;
		ChkPoint & operator=(const ChkPoint &) = delete;
		ChkPoint() = delete;
#endif

        struct Timestamp {
            const char * func;
            std::int32_t line;
            steady_clock::time_point realtime;
        };

        struct ChkPointFastImpl {
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
            static constexpr std::size_t N = 30;
#else
            static const std::size_t N = 30;
#endif
            std::int32_t cur;
            ChkPoint::Timestamp points[N];

            ChkPointFastImpl() : cur(0) {}
            ~ChkPointFastImpl() {}
        };

		template <typename T>
        struct fastpimpl_deleter {
            void operator()(T * const p) const {
                FastArenaObject<sizeof(ChkPoint::ChkPointFastImpl)>::
                    operator delete(reinterpret_cast<void *>(p));
            }
        };

		const std::unique_ptr<ChkPointFastImpl, fastpimpl_deleter<ChkPointFastImpl>> cfp;
	
	public:
#if (_MSC_VER >= 1700) || defined(__GXX_EXPERIMENTAL_CXX0X__)
		typedef std::pair<std::int64_t, std::int64_t> dpair;
#else
		typedef std::pair<const double, const double> dpair;
#endif
		ChkPoint(const char * const func, std::int32_t line);
		~ChkPoint();
		void checkpoint(const char * const func, std::int32_t line);
        ChkPoint::dpair totalpassageoftime() const;
		void checkpoint_print() const;
	};
}

#endif	// _CHKPOINT_H_
