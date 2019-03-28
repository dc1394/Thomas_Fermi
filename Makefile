#
# プログラム名
#
PROG = thomasfermi

#
# ソースコードが存在する相対パス
#
VPATH = src/alglib src/checkpoint src/thomasfermi src/thomasfermi/gausslegendre \
		src/thomasfermi/makerhoen src/thomasfermi/mixing src/thomasfermi/myfunctional \
		src/thomasfermi/shoot

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
CXXFLAGS = -Wall -Wextra -std=c++17 -xHOST -O3 -ipo -no-prec-div -pipe -I${MKLROOT}/include 

#
# リンク対象に含めるライブラリの指定
#
LDFLAGS = -L/home/dc1394/oss/boost_1_69_0/stage/icc/lib \
          -lboost_program_options \
		  -L${MKLROOT}/lib/intel64 \
		  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -ldl
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
