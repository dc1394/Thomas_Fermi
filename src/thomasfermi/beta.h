/*! \file beta.h
    \brief β(x)を計算するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _BETA_H_
#define _BETA_H_

#pragma once

#include <cstdint>	// for std::uint32_t
#include <vector>	// for std::vector

namespace thomasfermi {
	namespace femall {
        //! A class.
        /*!
            β(x)を計算するクラス
        */
		class Beta final
		{
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

            // #region メンバ関数

            //! A public member function (const).
            /*!
                operator()(double x)の宣言
                β(x)を計算して返す
                \param x xの値
                \return β(x)の値
            */
            double operator()(double x) const;

            // #endregion メンバ関数

        private:
            // #region メンバ変数

            //!  A private variable (constant).
            /*!
                メッシュが格納された動的配列のサイズ
            */
            std::size_t const size_;
            
            //!  A private variable (constant).
            /*!
                x方向のメッシュが格納された動的配列
            */
            std::vector<double> const xvec_;

            //!  A private variable (constant).
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

		inline double Beta::operator()(double x) const
		{
			auto klo = 0U;
			auto khi = static_cast<std::uint32_t>(size_ - 1);

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

			// yvec[i] = f(xvec[i]), yvec[i + 1] = f(xvec[i + 1])の二点を通る直線を代入
			return (yvec_[khi] - yvec_[klo]) / (xvec_[khi] - xvec_[klo]) * (x - xvec_[klo]) + yvec_[klo];
		}

		// #endregion メンバ関数
	}
}

#endif	// _BETA_H_
