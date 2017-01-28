#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = .
VPATH= $(HOME)/../RectArray $(HOME)/../ImageReconstruct $(HOME)/../Interpolator $(HOME)/../Filter 

CFLAGS = -ggdb 
# CFLAGS = -O3
CFLAGS += -I. -I$(HOME)/../RectArray/ -I$(HOME)/../ImageReconstruct -I$(HOME)/../Interpolator -I$(HOME)/../Filter
# INCLUDE = -I./../../Error

CXX =g++ -std=c++0x

SRCFILES:= $(wildcard $(HOME)/../RectArray/*.cpp $(HOME)/../ImageReconstruct/*.cpp $(HOME)/../Filter/*.cpp $(HOME)/../Interpolator/*.cpp )

OBJS:=$(patsubst %.cpp, %.o, $(SRCFILES))


%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

mainTest: GNUmakefile main.cpp ImageTest.cpp $(SRCFILES) $(OBJS)
	$(CXX) main.cpp ImageTest.cpp $(CFLAGS) $(OBJS) -o maintest.exe $(INCLUDE) -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_objdetect

clean:
	-@rm -rf *.o *.exe *.d

-include $(OBJS:.o=.d)
