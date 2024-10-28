SHELL=/bin/sh

PROGRAM=Main

CXX=g++
LD=g++

CFLAGS=-DDEBUG -g -ggdb
CXXFLAGS=$(CFLAGS) `soqt-config --cppflags`
LDFLAGS=-g `soqt-config --cppflags` `soqt-config --ldflags` `soqt-config --libs`


all: $(PROGRAM)


.SUFFIXES: .cpp .o

$(PROGRAM): Main.o Mesh.o Painter.o
#	$(LD)  $(LDFLAGS) $^ -o $(PROGRAM)
	$(LD)  $(LDFLAGS) -o $(PROGRAM) $^ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<


clean:
	rm -f $(PROGRAM)
	rm -f *.o
	rm -f *~
