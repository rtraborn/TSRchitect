CC=gcc
CFLAGS=-Wall -DWITHOUT_CAIRO
LDFLAGS=-lgenometools -lbam -lz
ifdef lib
  LDFLAGS += -L$(lib)
endif
GTHDR=/usr/local/include/genometools
SAMHDR=/usr/local/src/NGS-DIR/SAMTOOLS/samtools/

cagep:	cagep.c CagepNodeVisitor.c CagepNodeVisitor.h
	$(CC) $(CFLAGS) -I $(GTHDR) -I $(SAMHDR) -o cagep cagep.c CagepNodeVisitor.c $(LDFLAGS)

clean:	
	rm cagep
