TARGET1    = nestdaq-raw2tree-sglch
TARGET2    = nestdaq-raw2tree-pla

CXX       = g++
CXXFLAGS  = $(shell root-config --cflags) -std=c++17 -Werror -Wall -O2
LFLAGS    = 
ROOTLIBS  = $(shell root-config --libs)

all:	$(TARGET1) $(TARGET2)

$(TARGET1): $(TARGET1).o
	$(CXX) $(LFLAGS) -o $@ $^ $(ROOTLIBS)
$(TARGET2): $(TARGET2).o
	$(CXX) $(LFLAGS) -o $@ $^ $(ROOTLIBS)
%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $<
.PHONY : clean
clean:
	rm -rf $(TARGET1) $(TARGET2) *.o *~
