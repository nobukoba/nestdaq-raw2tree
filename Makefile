TARGET    = nestdaq-raw2tree-sglch

CXX       = g++
CXXFLAGS  = $(shell root-config --cflags) -std=c++17 -Werror -Wall -O2
LFLAGS    = 
ROOTLIBS  = $(shell root-config --libs)

all:	$(TARGET)

$(TARGET): $(TARGET).o
	$(CXX) $(LFLAGS) -o $@ $^ $(ROOTLIBS)
%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $<
.PHONY : clean
clean:
	rm -rf $(TARGET) *.o *~
