#ifndef _MAKERHOENERGY_H_
#define _MAKERHOENERGY_H_

#include "../Beta.h"
#include <cstdint>          // for std::int32_t
#include <fstream>          // for std::ofstream
#include <memory>           // for std::shared_ptr
#include <iomanip>          // for std::setprecision
#include <tuple>            // for std::tuple
#include <utility>          // for std::get
#include <boost/cast.hpp>   // for boost::numeric_cast

namespace Thomas_Fermi {
	namespace MakeRhoEn {
		class MakeRhoEnergy final
		{
            typedef std::tuple<std::shared_ptr<const FEM_ALL::Beta>,
                               std::vector<double>> parameter_type;

            // #region メンバ変数

            //! A static const variable.
            /*!
            エネルギーを計算するときの定数の値
            */
            static const double ALPHA;

            //! A private const variable.
            /*!
            x方向のメッシュ
            */
            const std::vector<double> xvec_;

            //! A private const variable.
            /*!
            x方向のメッシュの刻み幅
            */
            const double dx_;

            //! A private variable.
            /*!
            規格化のための定数
            */
            double s_;

            //! A private const variable.
            /*!
            ファイル出力するときのループの最大数
            */
            const std::int32_t max_;

            //! A private variable.
            /*!
            ファイル出力用のストリーム
            */
            std::ofstream ofs_;

            //! A private const variable.
            /*!
            ファイル出力用のストリーム
            */
            const double Z_;

            //! A private const variable.
            /*!
            Betaクラスのオブジェクトへのスマートポインタ
            */
            const std::shared_ptr<const FEM_ALL::Beta> pbeta_;

            //! A private const variable.
            /*!
            x方向のメッシュが格納されたstd::vectorのサイズ
            */
            const std::size_t size_;

            // #endregion メンバ変数

            // #region コンストラクタ

            //! デフォルトコンストラクタ（禁止）
            /*!
            */
            MakeRhoEnergy() = delete;

            //! コピーコンストラクタ（禁止）
            /*!
            \param コピー元のオブジェクト
            */
			MakeRhoEnergy(const MakeRhoEnergy &) = delete;

        public:
            //!  A constructor.
            /*!
            \param pt パラメータ
            \param Z 原子番号
            */
            MakeRhoEnergy(const parameter_type & pt, double Z);

            // #endregion コンストラクタ

            // #region operator=()

        private:
            //! operator=()（禁止）
            /*!
            \param コピー元のオブジェクト（禁止）
            \return コピー元のオブジェクト
            */
            MakeRhoEnergy & operator=(const MakeRhoEnergy &) = delete;

            // #endregion operator=()

            // #region メンバ関数

            //! 原子のエネルギーを求める
            /*!
            \return 原子のエネルギーを求める
            */
            double makeEnergy() const;

            //! 関数ρ~(x)の値を返す
            /*!
            \param x xの値
            \return ρ~(x)の値
            */
            double rhoTilde(double x) const;

        public:
            //! 計算結果をファイルに出力する
            /*!
            \param n Gauss-Legendreの分点
            \param usecilk Intel® Cilk™ Plusを使うかどうか
            \param usesimd SIMDを使うかどうか
            */
            void saveresult(std::size_t n, bool usecilk, bool usesimd);

        private:
            //! 関数ρ(x)の値をファイルに書き込む
            /*!
            \param filename 書き込むファイル名
            */
            void saverho(const std::string & filename);

            //! 関数ρ~(x)の値をファイルに書き込む
            /*!
            \param コピー元のオブジェクト（禁止）
            \return コピー元のオブジェクト
            */
            void saverhoTilde(const std::string & filename);

            //! 関数y(x)の値をファイルに書き込む
            /*!
            \param filename 書き込むファイル名
            */
            void savey(const std::string & filename);

            //! 関数y(x)の値を返す
            /*!
            \param x xの値
            \return y(x)の値
            */
            double y(double x) const;

            // #endregion メンバ関数
		};

		inline MakeRhoEnergy::MakeRhoEnergy(const parameter_type & pt, double Z)
			:	xvec_(std::get<1>(pt)),
                dx_((xvec_[2] - xvec_[1]) * 2.0),
                max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] / MakeRhoEnergy::ALPHA / dx_)),
                pbeta_(std::get<0>(pt)),
				size_(xvec_.size()),
                Z_(Z)
		{
            // 例外指定
            ofs_.exceptions(std::ios::badbit | std::ios::failbit);
            ofs_.setf(std::ios::fixed, std::ios::floatfield);
            ofs_ << std::setprecision(15);
		}

        //! x^yを計算する
        /*!
        \param x xの値    
        \param y yの値
        \return x^yの値
        */
		constexpr double power(double x, std::uint32_t y);
	}
}

#endif  // _MAKERHOENERGY_H_