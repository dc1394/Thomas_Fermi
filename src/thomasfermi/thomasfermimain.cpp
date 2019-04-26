/*! \file thomasfermimain.cpp
    \brief Thomas-Fermi方程式を有限要素法で解くコード
    Copyright © 2015-2019 @dc1394 All Rights Reserved.

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

#include "../checkpoint/checkpoint.h"
#include "getcomlineoption.h"
#include "goexit.h"
#include "iteration.h"
#include "makerhoen/makerhoenergy.h"
#include <cstdlib>                      // for EXIT_FAILURE, EXIT_SUCCESS
#include <iostream>                     // for std::cerr

int main(int argc, char * argv[])
{
    checkpoint::CheckPoint cp;

    thomasfermi::GetComLineOption mg;
    switch (mg.getopt(argc, argv)) {
    case -1:
        thomasfermi::goexit();

        return EXIT_FAILURE;
        break;

    case 0:
        break;

    case 1:
        thomasfermi::goexit();

        return EXIT_SUCCESS;
        break;

    default:
        BOOST_ASSERT(!"何かがおかしい！");
        break;
    }

    cp.checkpoint("処理開始", __LINE__);
    try {
        thomasfermi::femall::Iteration iter(mg.getpairdata());

        cp.checkpoint("初期関数生成処理", __LINE__);

        iter.Iterationloop();

        cp.checkpoint("Iterationループ処理", __LINE__);
        thomasfermi::makerhoen::MakeRhoEnergy mre(iter.PData()->gauss_legendre_integ_norm_, iter.makeresult(), iter.PData()->Z_);
        mre.saveresult();

        cp.checkpoint("結果出力処理", __LINE__);
    } catch (std::bad_alloc const &) {
        std::cerr << "メモリ確保に失敗しました。強制終了します。" << std::endl;
        thomasfermi::goexit();

        return EXIT_FAILURE;
    } catch (std::logic_error const & e) {
        std::cerr << e.what() << std::endl;
        thomasfermi::goexit();

        return EXIT_FAILURE;
    } catch (std::runtime_error const & e) {
        std::cerr << e.what() << std::endl;
        thomasfermi::goexit();

        return EXIT_FAILURE;
    }

    std::cout << "正常終了しました。結果をファイルに出力しました。\n";

    cp.checkpoint_print();
    cp.totalpassageoftime();
    checkpoint::usedmem();

    thomasfermi::goexit();

    return EXIT_SUCCESS;
}

