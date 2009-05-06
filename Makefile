#   File:       Makefile
#   Target:     lake
#   Sources:    lake.c recipes.c hunt.c stats.c toolbox.c


SOURCES  = lake.c hunt.c stats.c toolbox.c recipes.c
OBJECTS  = lake.o hunt.o stats.o toolbox.o recipes.o
INCLUDES = recipes.h
TEST_FILES  = xxx.smr xxx.inp xxx.dsm xxx.try 
TEST_FILES += lake.log test.gnu
TEST_FILES += work.gnuplot
TEST_FILES += plot.ps

ARCHIVE_FILES  = Makefile
ARCHIVE_FILES += $(SOURCES)
ARCHIVE_FILES += $(INCLUDES)
ARCHIVE_FILES += $(TEST_FILES)
ARCHIVE_FILES += lake


CC = gcc
LD_LIBS = -lm

lake: $(OBJECTS) Makefile
	gcc -o lake $(OBJECTS) $(LD_LIBS)

test: lake
	./lake < xxx.inp > xxx.log

all: lake archive

clean:
	/bin/rm -f *.o core *%

tar archive: lake.tar.gz

lake.tar.gz: $(SOURCES) $(ARCHIVE_FILES)
	rm -f  lake.tar.gz
	tar czf lake.tar.gz $(ARCHIVE_FILES)
	# gzip   lake.tar
	# mv lake.tar.gz lake.tgz

