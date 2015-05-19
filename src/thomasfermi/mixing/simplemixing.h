/*! \file simplemixing.h
	\brief 一次混合を行うクラスの宣言

Copyright ©  2015 @dc1394 All Rights Reserved.
This software is released under the BSD-2 License.
*/

#ifndef _SIMPLEMIXING_H_
#define _SIMPLEMIXING_H_

#pragma once

#include "../data.h"
#include <vector>

namespace thomasfermi {
	namespace mixing {
		class SimpleMixing final {
			// #region コンストラクタ・デストラクタ

			//! A constructor.
			/*!
				唯一のコンストラクタ
				\param arg インプットファイル名とCilk Plusを使用するかどうかのstd::pair
			*/
			SimpleMixing();

			//! A destructor.
			/*!
				デストラクタ
			*/
			~Iteration();
	}
}

#endif	// _SIMPLEMIXING_H_