#http://Evgenii.Rudnyi.Ru

include ../../make.inc

CPPFLAGS = -s -I.. -I../../../include 
LIB = ../minit.o 

%.o : %.cpp ../*.h 
	$(CC) -c $(CPPFLAGS) $<

all: mat.exe minmain.exe

mat.exe : mat.o
	$(CXX) -s $< -o $@
	$@

minmain.exe : minmain.o
	$(CXX) -s $< $(LIB) -o $@
	$@ -idat

clean:
	-rm *.o 

distclean:
	-rm *.o *.exe
