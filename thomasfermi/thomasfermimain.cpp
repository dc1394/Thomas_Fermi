#include "iteration.h"
#include "makerhoen/makerhoenergy.h"
#include <iostream>
#include <conio.h>

int main()
{
	//CheckPoint::ChkPoint cp("処理開始", __LINE__);
	try {
		thomasfermi::femall::Iteration iter(0.4, 0.001, 1000, 1.0E-12, true, true, 1.0E-5, 50.0, 5.0);

		//cp.checkpoint("初期関数生成処理", __LINE__);

		iter.Iterationloop();

		//cp.checkpoint("Iterationループ処理", __LINE__);
		thomasfermi::makerhoen::MakeRhoEnergy mre(1000, iter.makeresult(), true, 1);
		//mre.saveresult();

		//cp.checkpoint("結果出力処理", __LINE__);
	} catch (std::bad_alloc const &) {
		std::cerr << "メモリ確保に失敗しました。強制終了します。" << std::endl;

		std::cout << "終了するには何かキーを押してください..." ;
		::_getch();

		return EXIT_FAILURE;
	} catch (std::logic_error const & e) {
		std::cerr << e.what() << std::endl;

		std::cout << "終了するには何かキーを押してください..." ;
		::_getch();

		return EXIT_FAILURE;
	} catch (std::runtime_error const & e) {
		std::cerr << e.what() << std::endl;

		std::cout << "終了するには何かキーを押してください..." ;
		::_getch();

		return EXIT_FAILURE;
	}

	std::cout << "正常終了しました。結果をファイルに出力しました。\n";
	
	//cp.checkpoint_print();
	//CheckPoint::usedmem();

	std::cout << "終了するには何かキーを押してください..." ;
	::_getch();

	return EXIT_SUCCESS;
}
