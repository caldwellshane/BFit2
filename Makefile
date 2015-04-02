ifeq ($(ROOTSYS),)
$(error ROOTSYS is not define, do ssetup root)
endif
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
#INCDIRS := -I/opt/scarlet-3.x/include -I/opt/ktool-2.0/include -I/home/scaldwell/code/include 
CXX = g++
CXXFLAGS = -c -g -O2 $(ROOTCFLAGS) -I/opt/scarlet-3.x/include -I/opt/ktool-2.0/include -I/home/scaldwell/code/include #$(INCDIRS)
LIBS = -Wl,-rpath,/opt/scarlet-3.x/lib:/opt/ktool-2.0/lib \
 -L/opt/scarlet-3.x/lib -lscarlet -L/opt/ktool-2.0/lib -lktool

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $< -o $@

#%.d: %.cxx
#	$(CXX) -MM $(INCDIRS) $< -o $@

cxxsrcs = $(wildcard *.cxx)

.PHONY: all clean

targets = BFit2

all: $(targets)

BFit2: BFit2.o CSVtoStruct.o BFit2Model.o BFit2Populations.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)
	
-include $(cxxsrcs:.cxx=.d)

clean:
	rm -f $(targets) *.o
	
