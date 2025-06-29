# modified from Makefile in HsmlAndProject

# turns on warnings for if values to interpolate are out of interpolation bounds
#OPTIONS += -DBOUNDSWARNINGS

EXEC  = gethalomap.so

OBJS   =  gethalomap.o

INCL = Makefile

TEST = map_test

# Wall: warning option, fPIC: suitable for dynamic linking, 
# mcmodel: machine-specific, something to do with code size (copied from HsmlAndProject make)
# -D<bla>: adds #Define <bla> = 1 to the file
# ld creates executable file or library out of object files, sometimes -shared needed
# ld: dynamic linker
# -fopenmp: openmp version
# DEBUG: report on each value's interpolation and intermediate steps
OMP += -fopenmp
# OPTIONS += -DDEBUG
CFLAGS =  -O3 $(OPTIONS) $(OMP) -fPIC -Wall -mcmodel=medium
CFLAGsTEST = -O3 $(OPTIONS) $(OMP) -Wall -mcmodel=medium
LNKCMD =  $(CC) -shared $(OMP)

LIBS   =  -lm -g

CC     =  gcc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) -o $(EXEC) $(LIBS)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) $(TEST)

# test version of the code: normal executable; run to test the interpolators
test: $(OBJS)
	$(CC) -O3 $(OPTIONS) $(OBJS) $(OMP) -o $(TEST) $(LIBS)
    
