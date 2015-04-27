#include "../checkpoint/checkpoint.h"
#include "getcomlineoption.h"
#include "goexit.h"
#include "iteration.h"
#include "makerhoen/makerhoenergy.h"
#include <cstdlib>                      // for EXIT_FAILURE, EXIT_SUCCESS

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
		thomasfermi::makerhoen::MakeRhoEnergy mre(iter.PData()->gauss_legendre_integ_, iter.makeresult(), iter.PData()->Z_);
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
	checkpoint::usedmem();

	thomasfermi::goexit();
	
	return EXIT_SUCCESS;
}
