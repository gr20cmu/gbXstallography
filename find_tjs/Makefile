OBJS = main.o sub.o sub2.o symmetry.o


#  OBJSPS = ttex.o stex.o invpfmap.squ.o invpfmap.tri.o
#  OBJSFTN = srandf.o
#  OBJSC = rand.o drand.o srandc.o
#  OBJSEXT = randf.o drandf.o
 
SRCS=${OBJS:.o=.for}
#  SRCSPS=${OBJSPS:.o=.f}
#  SRCSFTN=${OBJSFTN:.o=.f}
#  SRCSC=${OBJSC:.o=.c}
#  SRCSEXT=${OBJSEXT:.o=.f}

MISC = Makefile README common.fi input.dat
#  MORESOURCE = 

LIBS=
FLAGS= -O3 -c -ffixed-line-length-none
#FLAGS= -g -c -O3 -w -s -static -fomit-frame-pointer -Wall -mpentiumpro \
#                -march=pentiumpro -malign-functions=4 -funroll-loops \
#                -fexpensive-optimizations -malign-double  \
#                -mwide-multiply
LFLAGS=     
FTNCOMP = gfortran
CCOMP=cc

find_tjs:   $(OBJS)
	$(FTNCOMP) $(LFLAGS) $(OBJS) $(LIBS) -o find_tjs

clean:
	rm -f $(OBJS)

archive:
	tar -cf 5param.tar $(SRCS) $(MISC) $(MORESOURCE)

main.o: main.f common.fi
	$(FTNCOMP) $(FLAGS) main.f 

sub.o: sub.f common.fi
	$(FTNCOMP) $(FLAGS) sub.f

sub2.o: sub2.f common.fi
	$(FTNCOMP) $(FLAGS) sub2.f

symmetry.o: symmetry.f common.fi
	$(FTNCOMP) $(FLAGS) symmetry.f 
