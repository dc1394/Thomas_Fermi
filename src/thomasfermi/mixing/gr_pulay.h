/*! \file gr_pulay.h
    \brief guaranteed-reduction Pulay法でyの合成を行うクラスの宣言
           Ref:
           D.R.Bowler and M.J.Gillan, Chem.Phys.Lett. 325, 475 (2000)

    Copyright ©  2015 @dc1394 All Rights Reserved.
    (but this is originally adapted by T.Ozaki for GR_Pulay.c from ADPACK)
    This software is released under the BSD-2 License.
*/

#ifndef _GR_PULAY_H_
#define _GR_PULAY_H_

#pragma once

#include "simplemixing.h"

namespace thomasfermi {
    namespace mixing {
        class GR_Pulay : public SimpleMixing {
            // #region コンストラクタ・デストラクタ

        public:
            //! A constructor.
            /*!
                唯一のコンストラクタ
                \param pdata データオブジェクト
            */
            GR_Pulay(std::shared_ptr<Data> const & pdata);

            //! A destructor.
            /*!
                デフォルトデストラクタ
            */
            ~GR_Pulay() = default;

            // #endregion コンストラクタ・デストラクタ

            // #region publicメンバ関数

            //! A public member function.
            /*!
                guaranteed-reduction Pulay法によって、yの合成を行う関数
                \param scfiter SCFの回数
                \param x xのメッシュ
                \param y 新しいy
                \return 合成後のy
            */
            femall::FEM::dmklvector operator()(std::int32_t scfiter, femall::FEM::dvector const & x, femall::FEM::dmklvector const & y);

            // #endregion publicメンバ関数

            // #region privateメンバ関数

        private:
            //! A private member variable (constant expression).
            /*!
                SCF Mixing history
            */
            static auto constexpr NUM_MIXING_PDM = 5;

            //! A private member function.
            /*!
                残差ノルムを求める関数
                \return 残差ノルム
            */
            femall::FEM::dmklvector getry();

            //! A private member function.
            /*!
                残差ノルムを求める関数
                \param newy 新しいy
                \param oldy 古いy
                \return 残差ノルム
            */
            femall::FEM::dmklvector getry(femall::FEM::dmklvector const & newy, femall::FEM::dmklvector const & oldy);

            //! A private member function.
            /*!
                yと残差ノルムの履歴を求める関数
                \param y 新しいy
            */
            void setyryarray();

            // #endregion publicメンバ関数

            // #region メンバ変数

            //! A private member variable.
            /*!
                残差ノルム
            */
            femall::FEM::dmklvector py;

            //! A private member variable.
            /*!
                残差ノルムの履歴を集めたstd::vector
            */
            std::vector<femall::FEM::dmklvector> ryarray_;
            
            //! A private member variable.
            /*!
                yの履歴を集めたstd::vector
            */
            std::vector<femall::FEM::dmklvector> yarray_;

            // #endregion メンバ変数

            // #region 禁止されたコンストラクタ・メンバ関数

            //! A private constructor (deleted).
            /*!
                デフォルトコンストラクタ（禁止）
            */
            GR_Pulay() = delete;

            //! A private copy constructor (deleted).
            /*!
                コピーコンストラクタ（禁止）
            */
            GR_Pulay(const GR_Pulay &) = delete;

            //! A private member function (deleted).
            /*!
                operator=()の宣言（禁止）
                \param コピー元のオブジェクト（未使用）
                \return コピー元のオブジェクト
            */
            GR_Pulay & operator=(const GR_Pulay &) = delete;

            // #endregion 禁止されたコンストラクタ・メンバ関数
        };
    }
}

#endif // _GR_PULAY_H_
