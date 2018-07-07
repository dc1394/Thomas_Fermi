/*! \file data.h
    \brief インプットファイルの各種データの構造体の宣言
    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _DATA_H_
#define _DATA_H_

#pragma once

#include "ci_string.h"
#include <array>                // for std::array
#include <cstdint>              // for std::uint32_t, std::uint8_t

namespace thomasfermi {
    //! A global variable (constant expression).
    /*!
        微分方程式を解くときの許容誤差のデフォルト値
    */
    static auto constexpr EPS_DEFAULT = 1.0E-15;

    //! A global variable (constant expression).
    /*!
        Gauss-Legendre積分の分点のデフォルト値
    */
    static auto constexpr GAUSS_LEGENDRE_INTEG_DEFAULT = 5U;

    //! A global variable (constant expression).
    /*!
        エネルギーを求めるときのGauss-Legendre積分の分点のデフォルト値
    */
    static auto constexpr GAUSS_LEGENDRE_INTEG_NORM_DEFAULT = 1000U;
    
    //! A global variable (constant expression).
    /*!
        微分方程式を解くときのメッシュの数のデフォルト値
    */
    static auto constexpr GRID_NUM_DEFAULT = 20000U;
    
    //! A global variable (constant expression).
    /*!
        マッチングポイント（xmin〜xmaxまでの間でなければならない）のデフォルト値
    */
    static auto constexpr MATCH_POINT_DEFAULT = 7.5;
        
    //! A global variable (constant expression).
    /*!
        ITERATIONの収束判定条件の値のデフォルト値
    */
    static auto constexpr ITERATION_CRITERION_DEFAULT = 1.0E-13;

    //! A global variable (constant expression).
    /*!
        ITERATIONの最大ループ回数のデフォルト値
    */
    static auto constexpr ITERATION_MAXITER_DEFAULT = 1000U;

    //! A global variable (constant expression).
    /*!
        電子密度を合成するときの重みのデフォルト値
    */
    static auto constexpr ITERATION_MIXING_WEIGHT_DEFAULT = 0.13;

    //! A global variable (constant expression).
    /*!
        微分方程式を解くときのメッシュの最大値のデフォルト値
    */
    static auto constexpr XMAX_DEFAULT = 100.0;

    //! A global variable (constant expression).
    /*!
        微分方程式を解くときのメッシュの最小値のデフォルト値
    */
    static auto constexpr XMIN_DEFAULT = 1.0E-5;

    //! A struct.
    /*!
        インプットファイルの各種データの構造体
    */
    struct Data final {
        // #region メンバ変数
        
        //!  A public member variable.
        /*!
            微分方程式を解くときの許容誤差
        */
        double eps_ = EPS_DEFAULT;

        //!  A public member variable.
        /*!
            Gauss-Legendre積分の分点
        */
        std::uint32_t gauss_legendre_integ_ = GAUSS_LEGENDRE_INTEG_DEFAULT;
        
        //!  A public member variable.
        /*!
            正規化を行う時のGauss-Legendre積分の分点
        */
        std::uint32_t gauss_legendre_integ_norm_ = GAUSS_LEGENDRE_INTEG_NORM_DEFAULT;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの数
        */
        std::uint32_t grid_num_ = GRID_NUM_DEFAULT;

        //!  A public member variable.
        /*!
            マッチングポイント（xmin〜xmaxまでの比率で表す）
        */
        double match_point_ = MATCH_POINT_DEFAULT;

        //!  A public member variable.
        /*!
            ITERATIONの収束判定条件の値
        */
        double iteration_criterion_ = ITERATION_CRITERION_DEFAULT;

        //!  A public member variable.
        /*!
            ITERATIONの最大ループ回数
        */
        std::uint32_t iteration_maxiter_ = ITERATION_MAXITER_DEFAULT;

        //!  A public member variable.
        /*!
            電子密度を合成するときの重み
        */
        double iteration_mixing_weight_ = ITERATION_MIXING_WEIGHT_DEFAULT;

        //!  A public member variable.
        /*!
            Cilk Plusを使用するかどうか
        */
        bool usecilk_ = true;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの最大値
        */
        double xmax_ = XMAX_DEFAULT;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの最小値
        */
        double xmin_ = XMIN_DEFAULT;

        //!  A public member variable.
        /*!
            原子核の電荷
        */
        double Z_;

        // #endregion メンバ変数
    };
}

#endif  // _DATA_H_

