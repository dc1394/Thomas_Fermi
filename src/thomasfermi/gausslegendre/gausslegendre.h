/*! \file Gauss_Legendre.h
    \brief Gauss-Legendre積分を行うクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD-2 License.
*/

#ifndef _GAUSS_LEGENDRE_H_
#define _GAUSS_LEGENDRE_H_

#pragma once

#include "../myfunctional/functional.h"
#include "../utility/property.h"
#include <array>                            // for std::array
#include <cstdint>                          // for std::int32_t
#include <vector>                           // for std::vector

namespace gausslegendre {
    //! A class.
    /*!
        Gauss-Legendre積分を行うクラス
    */
	class Gauss_Legendre final {
    public:
        // #region コンストラクタ・デストラクタ

        //! A constructor.
        /*!
            唯一のコンストラクタ
            Gauss-Legendreの重みと節を計算して、それぞれw_とx_に格納する
            \param n Gauss-Legendreの分点
        */
        explicit Gauss_Legendre(std::int32_t n);

		//! A destructor.
		/*!
			デフォルトデストラクタ
		*/
		~Gauss_Legendre() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

		template <typename FUNCTYPE>
        //! A public member function (template function).
        /*!
            Gauss-Legendre積分を実行する
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \return 積分値
        */
        double qgauss(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2) const;

    private:
        //! A private member function.
        /*!
            AVX命令が使用可能かどうかをチェックする
            \return AVX命令が使用可能ならtrue、使用不可能ならfalse
        */
        bool availableAVX() const;

        // #endregion メンバ関数

		// #region プロパティ

	public:
		//! A property.
		/*!

		*/
		utility::Property<std::vector<double>> const W;

		// #endregion プロパティ

        // #region メンバ変数

        //! A private member variable (constant).
        /*!
            Gauss-Legendreの分点
        */
        std::uint32_t const n_;

        //! A private member variable.
        /*!
			Gauss-Legendreの重み
        */
        std::vector<double> w_;

        //! A private member variable.
        /*!
            Gauss-Legendreの節（alignmentが揃っている）
        */
        std::vector<double> x_;
        
        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Gauss_Legendre() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
		Gauss_Legendre(Gauss_Legendre const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
		Gauss_Legendre & operator=(Gauss_Legendre const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

    // #region template publicメンバ関数

    template <typename FUNCTYPE>
    double Gauss_Legendre::qgauss(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2) const
    {
        auto const xm = 0.5 * (x1 + x2);
        auto const xr = 0.5 * (x2 - x1);

        auto sum = 0.0;
        for (auto i = 0U; i < n_; i++) {
			sum += w_[i] * func(xm + xr * x_[i]);
        }

        return sum * xr;
    }

    // #endregion template publicメンバ関数
}

#endif  // _GAUSS_LEGENDRE_H_