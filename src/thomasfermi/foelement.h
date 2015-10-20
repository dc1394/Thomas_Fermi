/*! \file foelement.h
    \brief 一次要素のクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _FOELEMENT_H_
#define _FOELEMENT_H_

#pragma once

#include "fem.h"

namespace thomasfermi {
    namespace femall {
        //! A class.
        /*!
            一次要素のクラス
        */
        class FOElement final : public FEM {
            // #region コンストラクタ・デストラクタ

        public:
            //! A constructor.
            /*!
                唯一のコンストラクタ
                \param beta
                \param coords
                \param nint
                \param usecilk Cilkを使用するかどうか
            */
            FOElement(dvector && beta, dvector const & coords, std::size_t nint, bool usecilk);

            //! A destructor.
            /*!
                デストラクタ
            */
            ~FOElement() override = default;

            // #endregion コンストラクタ・デストラクタ 

            // #region publicメンバ関数
            
            //! A public member function (constant - override).
            /*!
                結果を返す関数
                \return 結果を集めたboost::container::flat_map
            */
            FEM::resulttuple createresult() const override;

            //! A public member function (override).
            /*!
                βの状態をリセットする
                \param beta 対象のβ
            */
            void reset(dvector const & beta) override;

            // #endregion publicメンバ関数

            // #region privateメンバ関数

        private:
            //! A private member function.
            /*!
                a0_とa1_にastiff_を足し込む関数
                \param ielem
            */
            void amerge(std::size_t ielem) override;

            //! A private member function (pure virtual function).
            /*!
                小行列の要素を生成する
                \param ielem
            */
            void element(std::size_t ielem) override;
            
            //! A private member function (constant).
            /*!
                dn/drを返す関数
                \return dn/dr
            */
            dvector getdndr() const;

            //! A private member function (constant - override).
            /*!
                cを返す関数
                \param ielem 要素
                \return c
            */
            dvector getc(std::size_t ielem) const override;

            // #endregion メンバ関数

            // #region メンバ変数
            
            //! A private member variable.
            /*!
                double func(double)の形の関数オブジェクト
            */
            std::function<double(double)> func_;

            //! A private member variable (constant).
            /*!
                関数オブジェクト
            */
            std::function<double(double, double, std::size_t)> const fun1_;

            //! A private member variable (constant).
            /*!
                関数オブジェクト
            */
            std::function<double(double, double, std::size_t)> const fun2_;

            //! A private member variable.
            /*!
                形状関数1を格納する関数オブジェクト
            */
            std::function<double (double)> N1_;
            
            //! A private member variable.
            /*!
                形状関数2を格納する関数オブジェクト
            */
            std::function<double (double)> N2_;

            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            FOElement() = delete;

            //! A private copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
                \param コピー元のオブジェクト（未使用）
            */
            FOElement(FOElement const &) = delete;
            
            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            FOElement & operator=(FOElement const &) = delete;
            
            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }
}

#endif  // _FOELEMENT_H_
