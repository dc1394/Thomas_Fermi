#include "ChkPoint.h"
#include <iostream>
#include <boost/format.hpp>
#include <boost/assert.hpp>
#include <boost/optional.hpp>

#if (_MSC_VER < 1700) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
    #include <boost/chrono/chrono_io.hpp>
    #include <boost/chrono/process_cpu_clocks.hpp>
#endif

#ifdef _MSC_VER
	#include <windows.h>
	#include <Psapi.h>
#if (_MSC_VER >= 1700)
	#include <system_error>
#else
	#include <boost/system/error_code.hpp>
	#include <boost/system/system_error.hpp>
#endif
    #pragma comment(lib, "Psapi.lib")
#endif

namespace CheckPoint {
#if (_MSC_VER >= 1700) || defined(__GXX_EXPERIMENTAL_CXX0X__)
	using std::milli;
#else
	using boost::milli;
#endif

	ChkPoint::ChkPoint(const char * const func, std::int32_t line)
	 :	cfp(reinterpret_cast<ChkPointFastImpl *>(FastArenaObject<sizeof(ChkPointFastImpl)>::
												 operator new(0)))
	{
		checkpoint(func, line);
	}

	void ChkPoint::checkpoint(const char * const func, std::int32_t line)
	{
		BOOST_ASSERT(cfp->cur < ChkPoint::ChkPointFastImpl::N);

		Timestamp * const p = cfp->points + cfp->cur;
		p->func = func;
		p->line = line;
		p->realtime = steady_clock::now();

		cfp->cur++;
	}

	ChkPoint::dpair ChkPoint::totalpassageoftime() const
	{
		const auto realtime = duration_cast<milliseconds>(cfp->points[cfp->cur - 1].realtime -
			cfp->points[0].realtime);

		return std::make_pair(0, realtime.count());
	}

	void ChkPoint::checkpoint_print() const
	{
		boost::optional<steady_clock::time_point> prevreal(boost::none);

		for (std::int32_t i = 0; i < cfp->cur; i++) {
			Timestamp * const p = &cfp->points[i];

			if (prevreal) {
				const auto realtime(duration_cast<duration<double, milli>>(p->realtime - *prevreal));
				std::cout << p->func << " にかかった時間 = "
						  << boost::format("%.4f") % realtime.count()
						  << " (msec)\n";
			}

			prevreal = boost::optional<steady_clock::time_point>(p->realtime);
		}
	}

	ChkPoint::~ChkPoint()
	{
	}

#ifdef _MSC_VER
	void usedmem()
	{
		PROCESS_MEMORY_COUNTERS memInfo = {0};
		if (!::GetProcessMemoryInfo(::GetCurrentProcess(), &memInfo, sizeof(memInfo))) {
#if (_MSC_VER >= 1700)
			throw std::system_error(std::error_code(::GetLastError(), std::system_category()));
#else
			throw boost::system::system_error(boost::system::error_code(
											  ::GetLastError(),
											  boost::system::get_system_category()));
#endif
		}
		std::cout << "Used Memory Size: "
				  << static_cast<std::size_t>(memInfo.PeakWorkingSetSize >> 10)
				  << " (KiB)" << std::endl; 
	}
#endif
}
