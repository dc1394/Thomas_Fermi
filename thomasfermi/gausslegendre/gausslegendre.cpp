/*! \file Gauss_Legendre.cpp
    \brief Gauss-Legendre積分を行うクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/

#include "gausslegendre.h"
#include <stdexcept>        // for std::runtime_error
#include <boost/cast.hpp>   // for boost::numeric_cast

namespace gausslegendre {
    // #region コンストラクタ

    Gauss_Legendre::Gauss_Legendre(std::int32_t n)
        : avxSupported_(availableAVX()), n_(n)
    {
        alglib::ae_int_t info = 0;
        alglib::real_1d_array x, w;

        alglib::gqgenerategausslegendre(
            boost::numeric_cast<alglib::ae_int_t>(n),
            info,
            x,
            w);

        switch (info) {
        case 1:
            break;

        default:
            throw std::runtime_error("alglib::gqgenerategausslegendreが失敗");
            break;
        }

        x_.assign(x.getcontent(), x.getcontent() + x.length());
        w_.assign(w.getcontent(), w.getcontent() + w.length());
    }

    // #endregion コンストラクタ

    // #region privateメンバ関数

    bool Gauss_Legendre::availableAVX() const
    {
#if (_MSC_FULL_VER >= 160040219)
        std::array<std::int32_t, 4> cpuInfo = { 0 };
        ::__cpuid(cpuInfo.data(), 1);

        auto const osUsesXSAVE_XRSTORE = cpuInfo[2] & (1 << 27) || false;
        auto const cpuAVXSuport = cpuInfo[2] & (1 << 28) || false;

        if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
        {
            // Check if the OS will save the YMM registers
            auto const xcrFeatureMask = ::_xgetbv(_XCR_XFEATURE_ENABLED_MASK);
            return (xcrFeatureMask & 0x6) || false;
        }
#endif
        return false;
    }

    // #endregion privateメンバ関数
}
