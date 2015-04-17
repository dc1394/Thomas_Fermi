/*! \file MakeRhoEnergy.h
    \brief β(x)から電子密度とエネルギーを計算してファイルに記録するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _MAKERHOENERGY_H_
#define _MAKERHOENERGY_H_

#pragma once

#include "../Beta.h"
#include "../gausslegendre/gausslegendre.h"
#include <cstdint>                              // for std::int32_t
#include <fstream>                              // for std::ofstream
#include <memory>                               // for std::shared_ptr
#include <utility>                              // for std::pair

namespace thomasfermi {
	namespace makerhoen {
        //! A class.
        /*!
            y(x)から電子密度とエネルギーを計算する
        */
		class MakeRhoEnergy final {
            // #region 型エイリアス

            using parameter_type = std::pair < std::shared_ptr<fem_all::Beta>, std::vector<double> > ;

            // #endregion 型エイリアス

        public:
            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
				唯一のコンストラクタ
                \param n Gauss-Legendreの分点
                \param pt std::shared_ptr<Beta>とstd::vector<double>のタプル
                \param usesimd SIMDを使うかどうか
                \param Z 原子番号
            */
            MakeRhoEnergy(std::int32_t n, parameter_type const & pt, bool usesimd, double Z);

            //! A destructor.
            /*!
				デフォルトデストラクタ
            */
			~MakeRhoEnergy() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region publicメンバ関数

            //! A public member function.
            /*!
				計算結果をファイルに出力する
            */
            void saveresult();

            // #endregion publicメンバ関数

        private:
            // #region privateメンバ関数

            //! A private member function (const).
            /*!
                原子のエネルギーを求める
                \return 原子のエネルギーを求める
            */
            double makeEnergy() const;

            //! A private member function (const).
            /*!
                xを引数にとり、関数ρ(x)の値を返す
                \param x xの値
                \return ρ(x)の値
            */
            double rho(double x) const;

            //! A private member function (const).
            /*!
                xを引数にとり、関数ρ~(x)の値を返す
                \param x xの値
                \return ρ~(x)の値
            */
            double rhoTilde(double x) const;

            //! A private member function.
            /*!
                関数ρ(x)の値をファイルに書き込む
                \param filename 書き込むファイル名
            */
            void saverho(std::string const & filename);

            //! A private member function.
            /*!
                関数ρ~(x)の値をファイルに書き込む
                \param コピー元のオブジェクト（禁止）
                \return コピー元のオブジェクト
            */
            void saverhoTilde(std::string const & filename);

            //! A private member function.
            /*!
                関数y(x)の値をファイルに書き込む
                \param filename 書き込むファイル名
            */
            void savey(std::string const & filename);

            //! A private member function.
            /*!
                関数y(x)の値を返す
                \param x xの値
                \return y(x)の値
            */
            double y(double x) const;

            // #endregion privateメンバ関数

            // #region メンバ変数

            //! A private variable (constant).
            /*!
                エネルギーを計算するときの定数の値
            */
            double const alpha_;
			
			//! A private variable (constant).
			/*!
				x方向のメッシュが格納された動的配列
			*/
			std::vector<double> const xvec_;

			//! A private variable (constant).
			/*!
				x方向のメッシュの刻み幅
			*/
			double const dx_;

			//! A private variable (constant).
			/*!
				Gauss-Legendreの積分を行うオブジェクト
			*/
			gausslegendre::Gauss_Legendre const gl_;

			//! A private variable (constant).
			/*!
				ファイル出力するときのループの最大数
			*/
			std::int32_t const max_;

			//! A private variable (constant).
			/*!
				ファイル出力ストリーム
			*/
			std::ofstream ofs_;

			//! A private variable (constant).
			/*!
				Betaクラスのオブジェクトへのスマートポインタ
			*/
			std::shared_ptr<fem_all::Beta> const pbeta_;

			//! A private variable.
			/*!
				規格化のための定数
				s_ = 1.0 / (∫(0～∞)√x[y(x)]^(3/2)dx)
			*/
			double s_;

			//! A private variable (constant).
			/*!
				x方向のメッシュが格納された動的配列のサイズ
			*/
			std::size_t const size_;
			
			//! A private variable (constant).
			/*!
				SIMDを使うかどうか
			*/
			bool const usesimd_;
            
            //! A private variable (constant).
            /*!
                原子番号
            */
            double const Z_;
			            
            // #endregion メンバ変数

        private:
            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private default constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            MakeRhoEnergy() = delete;

            //! A private copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param コピー元のオブジェクト
            */
            MakeRhoEnergy(MakeRhoEnergy const &) = delete;

            //! operator=() (deleted).
            /*!
                operator=()（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            MakeRhoEnergy & operator=(MakeRhoEnergy const &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif  // _MAKERHOENERGY_H_
