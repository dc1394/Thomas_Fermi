PROG := thomasfermi
SRCS :=	alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp \
		fasttransforms.cpp integration.cpp interpolation.cpp linalg.cpp optimization.cpp \
		solvers.cpp specialfunctions.cpp specialfunctions.cpp statistics.cpp \
		checkpoint.cpp \
		beta.cpp ci_string.cpp fem.cpp foelement.cpp gausslegendre.cpp getcomlineoption.cpp \
		goexit.cpp iteration.cpp linearequations.cpp makerhoenergy.cpp readinputfile.cpp \
		simplemixing.cpp load2.cpp shootf.cpp shootfunc.cpp soelement.cpp thomasfermimain.cpp

OBJS :=	$(SRCS:%.cpp=%.o)
DEPS :=	$(SRCS:%.cpp=%.d)

VPATH  = src/alglib src/checkpoint src/thomasfermi src/thomasfermi/gausslegendre \
		 src/thomasfermi/makerhoen src/thomasfermi/mixing src/thomasfermi/shoot 
CXX = icpc
CXXFLAGS = -Wextra -O3 -pipe -std=c++14 -openmp -I${MKLROOT}/include
LDFLAGS = -L/home/dc1394/oss/boost_1_59_0/stage/icc/lib/ -lboost_program_options \
		  -lgsl -lgslcblas -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core \
		  -lmkl_intel_thread -lpthread

all: $(PROG) ;
#rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
