#http://Evgenii.Rudnyi.Ru

include ../make.inc

OBJ = enclofx.o \
      brent1.o  \
      chabis.o  \
      lbfgsb.o \
      tensolve.o \
      lbfgsbss.o

%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $<

../../include/%.h : %.h
	cp $? $@

../../bin/libtoms.a: $(OBJ) ../minit/minit.o
	$(ARCH) crs $@ $?

include: ../../include/toms.h

enclofx.o : toms.h
brent1.o : toms.h
chabis.o : toms.h
lbfgsbss.o : toms.h

lbfgsb.o : lbfgsb.f
	$(F77) $(F77FLAGS) -c $< 

tensolve.o : tensolve.f
	$(F77) $(F77FLAGS) -c $< 

clean:
	-rm *.o 

distclean:
	-rm *.o
