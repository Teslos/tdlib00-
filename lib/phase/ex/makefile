#http://Evgenii.Rudnyi.Ru

include ../../make.inc

CPPFLAGS = -s -I.. 
LIB1 = -L../../../bin -lphase 
LIB2 = -L../../../bin -lphase $(LAPACK) $(G2C)
LIB3 = -L../../../bin -lphase -ltoms $(LAPACK) $(G2C)


%.o : %.cpp ../*.h 
	$(CC) -c $(CPPFLAGS) $<

all: general calc coef elem formula func_tp species func_tpx interact phase simple cuox_pl
	
general: checkinf.exe enum.exe input.exe input2.exe tellg.exe

checkinf.exe : checkinf.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@

enum.exe : enum.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@

input.exe : input.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ input.mod

input2.exe : input2.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ input.mod

tellg.exe : tellg.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ 

calc: clc.exe
	
clc.exe : clc.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ 5+10*3.5
	$@ "exp(5)+10*sin(3.5)-7"

coef: cf.exe calc_coef.exe 
	
cf.exe : cf.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ cf.mod

calc_coef.exe : calc_coef.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ calc_coef.mod

elem: elem.exe
	
elm.exe : elm.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@

formula: frml.exe frml1.exe ins.exe

frml.exe : frml.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@

frml1.exe : frml1.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ frml1.mod

ins.exe : ins.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ ins.mod

func_tp: fnc_tp.exe
	
fnc_tp.exe : fnc_tp.o  
	$(CXX) -s $< $(LIB1) -o $@
	$@ fnc_tp.mod
	$@ fnc_tp_calc.mod
	$@ fnc_tp_cmpd.mod

species: sps.exe sps1.exe

sps.exe : sps.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ sps.mod

sps1.exe : sps1.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ sps.mod

func_tpx: fnc_tpx.exe fnc_tpx1.exe

fnc_tpx.exe : fnc_tpx.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ fnc_tpx.mod

fnc_tpx1.exe : fnc_tpx1.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ fnc_tpx1.mod

interact:  intrct.exe intrct2.exe intrct3.exe

intrct.exe : intrct.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ intrct2.mod

intrct2.exe : intrct2.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ intrct22.mod

intrct3.exe : intrct3.o  
	$(CXX) -s $< $(LIB2) -o $@
	$@ intrct2.mod

phase: phs.exe der.exe

phs.exe : phs.o  
	$(CXX) -s $< $(LIB3) -o $@
	$@ phs2.mod
	$@ ass.mod
	$@ y247.mod
	$@ y247_simple.mod

der.exe : der.o  
	$(CXX) -s $< $(LIB3) -o $@
	$@ phs2.mod der_2.dat

simple: iter.exe smpl.exe

iter.exe : iter.o  
	$(CXX) -s $< $(LIB3) -o $@
	$@

smpl.exe : smpl.o  
	$(CXX) -s $< $(LIB3) -o $@
	$@

cuox_pl: y123.exe

y123.exe : y123.o  
	$(CXX) -s $< $(LIB3) -o $@
	$@ y123.mod 1 -T 0.5 300 1000 100

clean:
	-rm *.o 

distclean:
	-rm *.o *.exe

