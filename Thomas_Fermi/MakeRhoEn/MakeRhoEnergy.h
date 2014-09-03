#include "../Beta.h"
#include <cstdint>
#include <fstream>
#include <memory>
#include <tuple>
#include <utility>
#include <boost/cast.hpp>

namespace Thomas_Fermi {
	namespace MakeRhoEn {
		class MakeRhoEnergy final
		{
            typedef std::tuple<std::vector<double>, std::shared_ptr<const FEM_ALL::Beta>,
                double> parameter_type;

            //! A static const variable.
            /*!
            エネルギーを計算するときの定数の値
            */
            static const double ALPHA;

            //! A private const variable.
            /*!
            x方向のメッシュの刻み幅
            */
            const double dx_;

            //! A private variable.
            /*!
            規格化のための定数
            */
            double k_;

            //! A private const variable.
            /*!
            ファイル出力するときのループの最大数
            */
            const std::int32_t max_;

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

            //! A private const variable.
            /*!
            x方向のメッシュ
            */
            const std::vector<double> xvec_;

            //! A private variable.
            /*!
            ファイル出力用のストリーム
            */
            std::ofstream ofs_;

            //! デフォルトコンストラクタ（禁止）
            /*!
              デフォルトコンストラクタ（禁止）
            */
            MakeRhoEnergy() = delete;

            //! コピーコンストラクタ（禁止）
            /*!
              \param コピー元のオブジェクト
            */
			MakeRhoEnergy(const MakeRhoEnergy &) = delete;

            //! operator=()（禁止）
            /*!
              \param コピー元のオブジェクト（禁止）
              \return コピー元のオブジェクト
            */
            MakeRhoEnergy & operator=(const MakeRhoEnergy &) = delete;

            //! 原子のエネルギーを求める
            /*!
            \return 原子のエネルギーを求める
            */
            double makeEnergy() const;

            //! 関数ρ(x)の値を返す
            /*!
            \param x xの値
            \return ρ(x)の値
            */
            double rhofunc(double x) const;

            //! 関数y(x)の値をファイルに書き込む
            /*!
            \param filename 書き込むファイル名
            */
            void savey(const std::string & filename);

            //! 関数ρ(x)の値をファイルに書き込む
            /*!
            \param filename 書き込むファイル名
            */
            void saverho(const std::string & filename);

            //! operator=()（禁止）
            /*!
            \param コピー元のオブジェクト（禁止）
            \return コピー元のオブジェクト
            */
            void saveresult3();

            //! 関数y(x)の値を返す
            /*!
              \param x xの値
              \return y(x)の値
            */
            double y(double x) const;







		public:
			MakeRhoEnergy::MakeRhoEnergy(const parameter_type & pt);

            //! 計算結果をファイルに出力する
            /*!
            \param n Gauss-Legendreの分点
            \param usecilk Intel Cilk Plusを使うかどうか
            \param usesimd SIMDを使うかどうか
            \param Z 原子番号
            */
			void saveresult(std::size_t n, bool usecilk, bool usesimd, double Z);
		};

		inline MakeRhoEnergy::MakeRhoEnergy(const parameter_type & pt)
			:	xvec_(std::get<0>(pt)), pbeta_(std::get<1>(pt)),
				size_(xvec_.size()),
				dx_((xvec_[2] - xvec_[1]) * 2.0),
				max_(boost::numeric_cast<std::int32_t>(xvec_[size_ - 1] /
					 MakeRhoEnergy::ALPHA / dx_))
		{
            // 例外指定
            ofs_.exceptions(std::ios::badbit | std::ios::failbit);
            ofs_.setf(std::ios::fixed, std::ios::floatfield);
            ofs_ << std::setprecision(15);
		}

		constexpr double power(double x, std::uint32_t y);
	}
}
