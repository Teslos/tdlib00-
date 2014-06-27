MAKE = make

all : include libs

include:
	(cd callback; $(MAKE) include)
	(cd f2c; $(MAKE) include)
	(cd lapack; $(MAKE) include)
	(cd minit; $(MAKE) include)
	(cd phase; $(MAKE) include)
	(cd td_algo; $(MAKE) include)
	(cd toms; $(MAKE) include)
	(cd varcomp; $(MAKE) include)

libs:
	(cd minit; $(MAKE))
	(cd toms; $(MAKE))
	(cd phase; $(MAKE))
	(cd td_algo; $(MAKE))
	(cd varcomp; $(MAKE))

distclean:
	rm -f ../include/*
	rm -f ../bin/*
	(cd minit; $(MAKE) distclean)
	(cd toms; $(MAKE) distclean)
	(cd phase; $(MAKE) distclean)
	(cd td_algo; $(MAKE) distclean)
	(cd varcomp; $(MAKE) distclean)

