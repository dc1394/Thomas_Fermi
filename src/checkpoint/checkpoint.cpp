/*! \file checkpoint.cpp
    \brief 時間計測のためのクラスの実装
    Copyright © 2014-2019 @dc1394 All Rights Reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "checkpoint.h"
#include <iostream>             // for std::cout
#include <optional>             // for std::nullopt, std::optional
#include <system_error>         // for std::system_category
#include <boost/assert.hpp>        // for boost::assert
#include <boost/cast.hpp>       // for boost::numeric_cast
#include <boost/format.hpp>     // for boost::format

#ifdef _WIN32
    #include <Windows.h>        // for GetCurrentProcess
    #include <Psapi.h>          // for GetProcessMemoryInfo

    #pragma comment(lib, "Psapi.Lib")
#else
    #include <errno.h>          // for errno 
    #include <sys/resource.h>   // for getrusage
#endif

namespace checkpoint {
    CheckPoint::CheckPoint()
        : cfp(
            reinterpret_cast<CheckPoint::CheckPointFastImpl *>(
                FastArenaObject<sizeof(CheckPoint::CheckPointFastImpl)>::operator new(0)))
    {
    }

    void CheckPoint::checkpoint(char const * action, std::int32_t line)
    {
        BOOST_ASSERT(cfp->cur < static_cast<std::int32_t>(CheckPoint::CheckPointFastImpl::N));

        auto const p = cfp->points.begin() + cfp->cur;

        p->action = action;
        p->line = line;
        p->realtime = std::chrono::high_resolution_clock::now();

        cfp->cur++;
    }
    
    void CheckPoint::checkpoint_print() const
    {
        using namespace std::chrono;

        std::optional<high_resolution_clock::time_point> prevreal(std::nullopt);

        auto itr = cfp->points.begin();
        for (auto i = 0; i < cfp->cur; ++i, ++itr) {
            if (prevreal) {
                auto const realtime(duration_cast< duration<double, std::milli> >(itr->realtime - *prevreal));
                std::cout << itr->action
                          << boost::format(" elapsed time = %.4f (msec)\n") % realtime.count();
            }

            prevreal = std::make_optional(itr->realtime);
        }
    }

    void CheckPoint::totalpassageoftime() const
    {
        using namespace std::chrono;

        auto const realtime = duration_cast< duration<double, std::milli> >(
            cfp->points[cfp->cur - 1].realtime - cfp->points[0].realtime);

        std::cout << boost::format("Total elapsed time = %.4f (msec)") % realtime.count() << std::endl;
    }

    // #region 非メンバ関数

#ifdef _WIN32
    void usedmem()
    {
        PROCESS_MEMORY_COUNTERS memInfo = { 0 };
        
        if (!::GetProcessMemoryInfo(::GetCurrentProcess(), &memInfo, sizeof(memInfo))) {
            throw std::system_error(std::error_code(::GetLastError(), std::system_category()));
        }

        std::cout << "Used Memory Size: "
                  << boost::numeric_cast<std::uint32_t>(memInfo.PeakWorkingSetSize >> 10)
                  << "(kB)"
                  << std::endl; 
    }
#else
    void usedmem()
    {
        struct rusage r;

        if (getrusage(RUSAGE_SELF, &r)) {
            throw std::system_error(errno, std::system_category());
        }

        std::cout << "Used Memory Size: "
                  << boost::numeric_cast<std::uint32_t>(r.ru_maxrss)
                  << "(kB)"
                  << std::endl;
    }
#endif

    // #endregion 非メンバ関数
}
