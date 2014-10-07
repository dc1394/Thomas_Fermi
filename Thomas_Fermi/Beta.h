/*! \file Beta.h
    \brief β(x)を計算するクラスの宣言

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#ifndef _BETA_H_
#define _BETA_H_

#pragma once

#include "mkl_allocator.h"
#include <vector>

namespace thomasfermi {
	namespace FEM_ALL {
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
                \param xvec x方向のメッシュ
                \param yvec y方向のメッシュ
            */
            Beta(const std::vector<double> & xvec, const std::vector<double> & yvec)
                : size_(xvec.size()), xvec_(xvec), yvec_(yvec) {}

            //! A destructor.
            /*!
                何もしないデストラクタ
            */
            ~Beta() {}

            // #endregion コンストラクタ・デストラクタ

        public:
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
			Beta(const Beta &) = delete;

            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
			Beta & operator=(const Beta &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
		};
	}
}

#endif	// _BETA_H_
