#http://Evgenii.Rudnyi.Ru/

include ../make.inc

ADDFLAGS = -I../phase

%.o : %.cpp
	$(CC) $(CPPFLAGS) $(ADDFLAGS) -c $<

../../include/%.h : %.h
	cp $? $@

OBJ = \
	common.o \
	data.o \
	residual.o \
	sumsqr.o \
	sumsqrpr.o \
	sumsqrut.o \
	sumsqrin.o \
	static.o \
	global.o \
	post_an.o \
	opt.o \
	hooke.o


INCS = common.h \
	data.h \
	residual.h \
	sumsqr.h \
	post_an.h

all : ../../bin/libvarcomp.a ../../bin/assess

../../bin/libvarcomp.a: $(OBJ)
	$(ARCH) crs $@ $?

../../bin/assess: assess.o ../../bin/libvarcomp.a ../../bin/libtd_algo.a ../../bin/libphase.a ../../bin/libtoms.a
	$(CXX) -s $(LDFLAGS) $< -o $@ -L../../bin -lvarcomp -ltd_algo -lphase -ltoms $(LAPACK) $(ADDLIB) $(G2C)

# -O3 breaks the hooke
hooke.o : hooke.cpp
	$(CC) $(CPPFLAGS) $(ADDFLAGS) -O0 -c $<

common.o : common.h
data.o : data.h common.h
residual.o : sumsqr.h residual.h data.h common.h ../../include/td_algo.h
sumsqr.o : opt.h sumsqr.h residual.h data.h common.h ../../include/td_algo.h
sumsqrin.o : sumsqr.h residual.h data.h common.h ../../include/td_algo.h
sumsqrpr.o : sumsqr.h residual.h data.h common.h ../../include/td_algo.h
sumsqrut.o : sumsqr.h residual.h data.h common.h ../../include/td_algo.h
static.o : sumsqr.h residual.h data.h common.h ../../include/td_algo.h
global.o : sumsqr.h data.h common.h ../../include/td_algo.h
post_an.o : sumsqr.h post_an.h data.h common.h ../../include/td_algo.h
opt.o : opt.h
hooke.o : opt.h

assess.o : sumsqr.h residual.h data.h common.h ../../include/td_algo.h

include: $(addprefix ../../include/, $(INCS))

clean:
	-rm *.o

distclean:
	-rm *.o
