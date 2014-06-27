#http://Evgenii.Rudnyi.Ru

include ../make.inc

OBJ = minit.o 

%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $<

../../include/%.h : %.h
	cp $? $@

minit.o: minit.cpp

include: ../../include/minit.h

minit.o : minit.h

clean:
	-rm *.o 

distclean:
	-rm *.o 
