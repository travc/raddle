DEST = ../bin/
CC = gcc
INCDIRS = -I/opt/local/include
LIBDIRS = -L/opt/local/lib

OBJS =          ranlib.o  \
                com.o     \
                linpack.o\
                ECA_utilities.o  \
                ECA_MemAlloc.o  \
                MathStatRand.o  \
		JN_RVutils.o \
		JN_MemAlloc.o \
		JN_TreeUtils.o \
		JN_MCTypes.o \
		RaddleBase.o


GOPT = -O3
HFLAGS=
CFLAGS = $(GOPT) $(INCDIRS) $(HFLAGS) -Wall
LFLAGS = -lm -lgsl -lgslcblas $(LIBDIRS)
SRCS =  $(OBJS:.o=.c)


RaddleMC : $(OBJS) mainMC.o RaddleMCutils.o 
	$(CC) $(CFLAGS) -o $(DEST)/RaddleMC $(OBJS) mainMC.o RaddleMCutils.o $(LFLAGS) 

clean:
	rm -f *.o


#/usr/lib64/libgsl.a /usr/lib64/libgslcblas.a 
