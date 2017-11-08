#
# プログラム名
#
PROG = thomasfermi

#
# ソースコードが存在する相対パス
#
VPATH = src/alglib src/checkpoint src/thomasfermi src/thomasfermi/gausslegendre \
		src/thomasfermi/makerhoen src/thomasfermi/mixing src/thomasfermi/myfunctional \
		src/thomasfermi/shoot src/tho

#
# コンパイル対象のソースファイル群（カレントディレクトリ以下の*.cppファイル）
#
SRCS = $(shell find * -name "*.cpp")

#
# ターゲットファイルを生成するために利用するオブジェクトファイル
#
OBJDIR = 
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR = .
endif

OBJS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cpp=.o)))

#
# *.cppファイルの依存関係が書かれた*.dファイル
#
DEPS = $(OBJS:.o=.d)

#
# C++コンパイラの指定
#
CXX = icpc

#
# C++コンパイラに与える、（最適化等の）オプション
#
CXXFLAGS = -Wextra -O3 -xCORE-AVX512 -ip -no-prec-div -DMKL_ILP64 -pipe -std=c++17 -I/home/dc1394/openblas/include 

#
# リンク対象に含めるライブラリの指定
#
LDFLAGS = -L/home/dc1394/oss/boost_1_65_1/stage/icc/lib -L/home/dc1394/openblas/lib \
          -lboost_program_options -ldl -lgsl -lgslcblas -lifcore -liomp5 -lm -lopenblas \
		  -lpthread -lcilkrts
#
# makeの動作
#
all: $(PROG) ;

#
# 依存関係を解決するためのinclude文
#
-include $(DEPS)

#
# プログラムのリンク
#
$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

#
# プログラムのコンパイル
#
%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

#
# make cleanの動作
#
clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
