#include "checkpoint.h"
#include "goexit.h"
#include "iteration.h"
#include "makerhoen/makerhoenergy.h"
#include <cstdlib>                      // for EXIT_FAILURE, EXIT_SUCCESS

int main()
{
	checkpoint::CheckPoint cp;
	cp.checkpoint("処理開始", __LINE__);
	try {
		thomasfermi::femall::Iteration iter(1.0, 0.001, 1000, 1.0E-12, true, true, 1.0E-5, 50.0, 10.0);

		cp.checkpoint("初期関数生成処理", __LINE__);

		iter.Iterationloop();

		cp.checkpoint("Iterationループ処理", __LINE__);
        thomasfermi::makerhoen::MakeRhoEnergy mre(1000, iter.makeresult(), true, 1);
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
