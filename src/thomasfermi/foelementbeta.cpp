/*! \file foelementbeta.h
    \brief 一次要素でβ(x)を計算するクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "foelementbeta.h"

namespace thomasfermi {
    namespace femall {
        // #region コンストラクタ

        FoelementBeta::FoelementBeta(std::shared_ptr<Beta> const & pbeta) :
            pbeta_(pbeta)
        {
        }

        // #endregion コンストラクタ
    }
}
