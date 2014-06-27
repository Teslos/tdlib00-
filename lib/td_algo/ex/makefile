#http://Evgenii.Rudnyi.Ru/

include ../../make.inc

CPPFLAGS = -s -I.. -I../../../include -I/sw/include
LIB = -L../../../bin -ltd_algo -lvarcomp -lphase -ltoms $(LAPACK) $(G2C)

%.o : %.cpp ../*.h 
	$(CC) -c $(CPPFLAGS) $<

all: convert.exe evaluate.exe
	
convert.exe : convert.o  
	$(CXX) -s $< $(LIB) -o $@
	$@ cv.mod

#does not work properly: to fix
evaluate.exe : evaluate.o  
	$(CXX) -s $< $(LIB) -o $@
	-$@

clean:
	-rm *.o 

distclean:
	-rm *.o *.exe

