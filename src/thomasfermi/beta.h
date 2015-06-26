/*! \file beta.h
    \brief β(x)を計算するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _BETA_H_
#define _BETA_H_

#pragma once

#include "utility/property.h"
#include <cstdint>				// for std::uint32_t
#include <vector>				// for std::vector

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			β(x)を計算するクラス
		*/
		class Beta final {
        public:
            // #region コンストラクタ・デストラクタ

            //! A constructor.
            /*!
				コンストラクタ
                \param xvec x方向のメッシュ
                \param yvec y方向のメッシュ
            */
			Beta(std::vector<double> const & xvec, std::vector<double> const & yvec);
			
            //! A destructor.
            /*!
                デフォルトデストラクタ
            */
			~Beta() = default;

            // #endregion コンストラクタ・デストラクタ

			// #region プロパティ

			//! A property.
			/*!
				メッシュが格納された動的配列のサイズへのプロパティ
			*/
			utility::Property<std::size_t> const Size;
			
			//! A property.
			/*!
				x方向のメッシュが格納された動的配列へのプロパティ
			*/
			utility::Property<std::vector<double>> const Xvec;

			//! A property.
			/*!
				y方向のメッシュが格納された動的配列へのプロパティ
			*/
			utility::Property<std::vector<double>> const Yvec;

        private:
            // #region メンバ変数

            //!  A private member variable (constant).
            /*!
                メッシュが格納された動的配列のサイズ
            */
            std::size_t const size_;
            
            //!  A private member variable (constant).
            /*!
                x方向のメッシュが格納された動的配列
            */
            std::vector<double> const xvec_;

            //!  A private member variable (constant).
            /*!
				y方向のメッシュが格納された動的配列
            */
            std::vector<double> const yvec_;

            // #endregion メンバ変数

            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private constructor (deleted).
            /*!
				デフォルトコンストラクタ（禁止）
            */
            Beta() = delete;

            //! A private copy constructor (deleted).
            /*!
				コピーコンストラクタ（禁止）
            */
			Beta(Beta const &) = delete;

            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
			Beta & operator=(Beta const &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
		};

		// #region メンバ関数

		

        

        //    // yvec[i] = f(xvec[i]), yvec[i + 1] = f(xvec[i + 1]), yvec[i + 2] = f(xvec[i + 2])の三点を通る放物線を生成
        //    
        //    // もし、配列の外にはみ出るときは
        //    if (khi >= max - 1) {
        //        // 一つ前の値を使う
        //        khi--;
        //        klo--;
        //    }

        //    auto const x2mx1 = xvec_[khi + 1] - xvec_[khi];
        //    auto const x0mx2 = xvec_[klo] - xvec_[khi + 1];
        //    auto const x1mx0 = xvec_[khi] - xvec_[klo];

        //    auto const denom = x2mx1 * x0mx2 * x1mx0;

        //    auto const a = -(x2mx1 * yvec_[klo] + x0mx2 * yvec_[khi] + x1mx0 * yvec_[khi + 1]) / denom;
        //    auto b = x2mx1 * (xvec_[khi + 1] + xvec_[khi]) * yvec_[klo];
        //    b += x0mx2 * (xvec_[klo] + xvec_[khi + 1]) * yvec_[khi];
        //    b += x1mx0 * (xvec_[khi] + xvec_[klo]) * yvec_[khi + 1];
        //    b /= denom;


        //    

        //}
		// #endregion メンバ関数
	}
}

#endif	// _BETA_H_
