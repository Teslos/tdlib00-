#http://Evgenii.Rudnyi.Ru/

include ../../make.inc

CPPFLAGS = -s -I.. -I../../../include

%.o : %.cpp ../*.h 
	$(CC) -c $(CPPFLAGS) $<

all: tst1.exe invert1.exe invert2.exe gauss1.exe
	
tst1.exe : tst1.o 
	$(CXX) -s $< $(LAPACK) $(G2C) -o $@
	$@

invert1.exe : invert1.o 
	$(CXX) -s $< $(LAPACK) $(G2C) -o $@
	$@

invert2.exe : invert2.o 
	$(CXX) -s $< $(LAPACK) $(G2C) -o $@
	$@

gauss1.exe : gauss1.o 
	$(CXX) -s $< $(LAPACK) $(G2C) -o $@
	$@ 100

clean:
	-rm *.o 

distclean:
	-rm *.o *.exe
