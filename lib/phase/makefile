#http://Evgenii.Rudnyi.Ru/

include ../make.inc

%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $<

../../include/%.h : %.h
	cp $? $@

OBJ = general.o \
      calc.o    \
      coef.o    \
      elem.o    \
      formula.o \
      func_tp.o \
      species.o \
      func_tpx.o \
      interact.o \
      phase.o   \
      simple.o  \
      cuox_or.o \
      intvars.o \
      cuox_pl.o \
      ass_sol.o 

INCS = \
 calc.h \
 coef.h \
 elem.h \
 formula.h \
 func_tp.h \
 func_tpx.h \
 general.h \
 phase.h \
 species.h 

../../bin/libphase.a: $(OBJ)
	$(ARCH) crs $@ $?

general.o : general.h
calc.o: calc.h general.h
coef.o: coef.h calc.h general.h
elem.o: elem.h general.h
formula.o: formula.h elem.h general.h
func_tp.o: func_tp.h func_tp_imp.h coef.h calc.h general.h
species.o: species.h func_tp.h func_tp_imp.h formula.h elem.h coef.h calc.h general.h ../../include/lapack.h
func_tpx.o: func_tpx.h func_tpx_imp.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h
interact.o: interact.h func_tpx.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h 
phase.o: phase.h func_tpx.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h
intvars.o: intvars.h phase.h func_tpx.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h ../../include/lapack.h
simple.o: simple.h phase.h func_tpx.h func_tpx_imp.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h
cuox_or.o: cuox_or.h phase.h func_tpx.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h
cuox_pl.o: cuox_pl.h intvars.h phase.h func_tpx.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h
apbq.o: apbq.h intvars.h simple.h phase.h func_tpx.h func_tpx_imp.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h
ass_sol.o: ass_sol.h intvars.h simple.h phase.h func_tpx.h func_tpx_imp.h species.h func_tp.h formula.h elem.h coef.h calc.h general.h ../../include/lapack.h ../../include/toms.h

include: $(addprefix ../../include/, $(INCS))

clean:
	-rm *.o

distclean:
	-rm *.o
