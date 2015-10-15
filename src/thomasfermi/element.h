/*! \file element.h
    \brief 要素を表す列挙型の宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#pragma once

#include <cstdint>  // for std::int32_t

namespace thomasfermi {
    namespace femall {
        enum class Element : std::int32_t {
            First = 1,
            Second = 2
        };
    }
}

#endif  // _ELEMENT_H_
