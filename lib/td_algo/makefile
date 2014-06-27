#http://Evgenii.Rudnyi.Ru/

include ../make.inc

%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $<

../../include/%.h : %.h
	cp $? $@

OBJ = td_algo.o \
  output.o \
  phase_eq.o

INCS =  td_algo.h \
  output.h 

../../bin/libtd_algo.a: $(OBJ)
	$(ARCH) crs $@ $?

td_algo.o : td_algo.h td_algo_imp.h ../../include/general.h ../../include/phase.h ../../include/coef.h
output.o : output.h td_algo.h ../../include/general.h
phase_eq.o : phase_eq.h td_algo.h ../../include/general.h ../../include/toms.h
evaluate.o : output.h td_algo.h ../../include/general.h

include: $(addprefix ../../include/, $(INCS))

clean:
	-rm *.o

distclean:
	-rm *.o

