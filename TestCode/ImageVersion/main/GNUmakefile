#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = .
VPATH= $(HOME)/../RectArray $(HOME)/../Image $(HOME)/../ImageReconstruct $(HOME)/../Interpolator $(HOME)/../Filter 

#CFLAGS = -ggdb 
CFLAGS = -O3
CFLAGS += -I. -I$(HOME)/../RectArray/ -I$(HOME)/../ImageReconstruct -I$(HOME)/../Interpolator -I$(HOME)/../Filter -I$(HOME)/../Image
# INCLUDE = -I./../../Error

SUBDIRS =.  $(HOME)/../RectArray $(HOME)/../Image $(HOME)/../ImageReconstruct $(HOME)/../Interpolator $(HOME)/../Filter

CXX =g++ -std=c++0x

SRCFILES:= $(wildcard $(HOME)/../RectArray/*.cpp $(HOME)/../ImageReconstruct/*.cpp $(HOME)/../Filter/*.cpp $(HOME)/../Interpolator/*.cpp $(HOME)/../Image/*.cpp)

OBJS:=$(patsubst %.cpp, %.o, $(SRCFILES))


%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

mainTest: GNUmakefile main.cpp $(SRCFILES) $(OBJS)
	$(CXX) main.cpp $(CFLAGS) $(OBJS) -o maintest.exe 

clean:
	-@rm -rf *.o *.exe *.d test*.pgm
moreclean:
	@for d in $(SUBDIRS); do (echo "Clearing in $$d/... "; cd $$d; $(MAKE) -s clean); done
-include $(OBJS:.o=.d)
