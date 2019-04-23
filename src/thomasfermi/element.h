/*! \file element.h
    \brief 要素を表す列挙型の宣言
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
