﻿#include "Iteration.h"
#include "MakeRhoEn/MakeRhoEnergy.h"
#include "ChkPoint.h"
#include <iostream>
#include <conio.h>

int main()
{
	CheckPoint::ChkPoint cp("処理開始", __LINE__);
	try {
		Thomas_Fermi::FEM_ALL::Iteration scf(1.0E-5, 20.0, 5.0, 0.001, 1000, false, true, 1.0E-12, 0.3);

		cp.checkpoint("初期関数生成処理", __LINE__);

		scf.Iterationloop();

		cp.checkpoint("Iterationループ処理", __LINE__);

		Thomas_Fermi::MakeRhoEn::MakeRhoEnergy mre(scf.makeresult());
		mre.saveresult(1000, true, true);

		cp.checkpoint("結果出力処理", __LINE__);
	} catch (const std::bad_alloc &) {
		std::cerr << "メモリ確保に失敗しました。強制終了します。" << std::endl;

		std::cout << "終了するには何かキーを押してください..." ;
		::_getch();

		return EXIT_FAILURE;
	} catch (const std::logic_error & e) {
		std::cerr << e.what() << std::endl;

		std::cout << "終了するには何かキーを押してください..." ;
		::_getch();

		return EXIT_FAILURE;
	} catch (const std::runtime_error & e) {
		std::cerr << e.what() << std::endl;

		std::cout << "終了するには何かキーを押してください..." ;
		::_getch();

		return EXIT_FAILURE;
	}

	std::cout << "正常終了しました。結果をファイルに出力しました。\n";
	
	cp.checkpoint_print();
	CheckPoint::usedmem();

	std::cout << "終了するには何かキーを押してください..." ;
	::_getch();

	return EXIT_SUCCESS;
}