PROJECT=int2         # for historical reasons we always like to 
                     # call the executable int2, but you could set it to 
                     # something more descriptive


OBJSUF = o
FC = gfortran -g -w -c
FFLAGS = -O2
FLINK = gfortran -o $(PROJECT)



.PHONY: $(PROJECT) clean list pack

.f.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.SUFFIXES: .$(OBJSUF) .f .c

# SOURCE FILE LIST
#
vpath %.f .:../src:../contrib/utilities/src


FSRCS = emdyadic_dr.f emdyadic.f emplanew.f prini.f

ifeq ($(WITH_SECOND),1) 
FSRCS += second-r8.f
endif

#
# object files list
OBJS    =  $(FSRCS:.f=.$(OBJSUF)) 
#
$(PROJECT): $(OBJS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJS)
	./$(PROJECT)
#
clean: 
	rm -f $(OBJS)
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
