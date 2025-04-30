
GCC = g++
SRCDIR = src
BUILDDIR = .build
BINDIR = bin
INCDIR = include
TARGET = $(BINDIR)/rsRNASP1

SRCPDB = $(SRCDIR)/PDB.cc \
         $(SRCDIR)/pdb_utils.cc
OBJ_PDB = $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SRCPDB:.cc=.o))

SRCrsRNASP = $(SRCDIR)/rsRNASP_calculator.cc \
           $(SRCDIR)/rsRNASP_PDB.cc \
           $(SRCDIR)/rsRNASP_dihedral.cc
OBJ_rsRNASP = $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SRCrsRNASP:.cc=.o))

CFLAGS = -std=c++11 -fopenmp# -g -pg
LDFLAGS = -Llib -fopenmp

LIB = -ldl 

INC = -I include


$(TARGET): $(OBJ_PDB) $(OBJ_rsRNASP) $(SRCDIR)/main.cc $(INCDIR)/main.h
	@mkdir -p $(BINDIR)
	@echo " Linking..." 
	$(GCC) $^ -o $(TARGET) $(LIB) $(LDFLAGS) $(INC) $(CFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cc $(INCDIR)/%.h
	@mkdir -p $(BUILDDIR)
	$(GCC) $(CFLAGS) $(INC) -c -o $@ $<

.PHONY: clean
clean:
	@echo " Cleaning..."
	rm -rf $(BUILDDIR) $(BINDIR)

# For training
train: $(OBJ_PDB) $(OBJ_rsRNASP) src/train.cc
	@mkdir -p $(BINDIR)
	$(GCC) $^ $(CFLAGS) $(INC) -static $(LIB) $(LDFLAGS)  -o bin/train

# get fasta sequences from PDB
getseq: $(OBJ_PDB) src/getseq.cc
	@mkdir -p $(BINDIR)
	$(GCC) $^ $(CFLAGS)  $(INC)  $(LIB) $(LDFLAGS)  -o bin/getseq
