#http://Evgenii.Rudnyi.Ru

include ../../make.inc

CPPFLAGS = -s -I.. -I../../../include -I/sw/include
LIB = -L../../../bin -ltoms -lphase $(G2C)

%.o : %.cpp ../*.h 
	$(CC) -c $(CPPFLAGS) $<

all: t_brent1.exe t_brent1a.exe t_chabis.exe t_chabisa.exe t_lbfgsb.exe t_lbfgsba.exe t_lbfgsbb.exe

t_brent1.exe : t_brent1.o  
	$(CXX) -s $< $(LIB) -o $@
	$@

t_brent1a.exe : t_brent1a.o  
	$(CXX) -s $< $(LIB) -o $@
	$@

t_chabis.exe : t_chabis.o  
	$(CXX) -s $< $(LIB) $(LAPACK) -o $@
	$@

t_chabisa.exe : t_chabisa.o  
	$(CXX) -s $< $(LIB) $(LAPACK) -o $@
	$@

t_lbfgsb.exe : t_lbfgsb.o  
	$(CXX) -s $< $(LIB) $(LAPACK) -o $@
	$@

t_lbfgsba.exe : t_lbfgsba.o  
	$(CXX) -s $< $(LIB) $(LAPACK) -o $@
	$@

t_lbfgsbb.exe : t_lbfgsbb.o  
	$(CXX) -s $< $(LIB) $(LAPACK) -o $@
	$@

clean:
	-rm *.o 

distclean:
	-rm *.o *.exe debug*

