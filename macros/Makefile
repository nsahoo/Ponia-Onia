#GCCFLAGS=-Wall -g
#ROOTFLAGS='root-config --cflags --libs'


GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

#TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTFLAGS) -c 
#SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
#FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS) 
#FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

PLHISTFLAGS=$(GCCFLAGS) -I$(ROOTFLAGS) -c
EFFHISTFLAGS=$(GCCFLAGS) -I$(ROOTFLAGS) -c

all: plothisto eff9 
	@echo "Full build successful."

plothisto: plotHisto.cc
	g++ $(PLHISTFLAGS) $< -o $@

eff9: eff_9aug.cc
	g++ $(EFFHISTFLAGS) $< -o $@

clean: 
	rm -rf plothisto eff9

