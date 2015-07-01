/*! \file beta.h
    \brief β(x)を計算するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _BETA_H_
#define _BETA_H_

#pragma once

#include "element.h"
#include "utility/property.h"
#include <cstdint>				// for std::uint32_t
#include <memory>               // for std::shared_ptr
#include <vector>				// for std::vector

namespace thomasfermi {
	namespace femall {
		//! A class.
		/*!
			β(x)を計算するクラス
		*/
		class Beta final {
            // #region コンストラクタ・デストラクタ

        public:
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

            // #region メンバ関数

			template <Element E>
            //! A public member function (const).
            /*!
                operator()(double x)の宣言
                一次要素でβ(x)を計算して返す
                \param x xの値
                \return β(x)の値
            */
            double operator()(double x) const;

            // #endregion メンバ関数

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

		template <>
		double Beta::operator()<Element::First>(double x) const
		{
			auto klo = 0U;
			auto const max = static_cast<std::uint32_t>(size_ - 1);
			auto khi = max;

			// 表の中の正しい位置を二分探索で求める
			while (khi - klo > 1) {
				auto const k = static_cast<std::uint32_t>((khi + klo) >> 1);

				if (xvec_[k] > x) {
					khi = k;
				}
				else {
					klo = k;
				}
			}

			// yvec_[i] = f(xvec_[i]), yvec_[i + 1] = f(xvec_[i + 1])の二点を通る直線を代入
			return (yvec_[khi] - yvec_[klo]) / (xvec_[khi] - xvec_[klo]) * (x - xvec_[klo]) + yvec_[klo];
		}
		       

        //    // yvec_[i] = f(xvec_[i]), yvec_[i + 1] = f(xvec_[i + 1]), yvec_[i + 2] = f(xvec_[i + 2])の三点を通る放物線を生成
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
