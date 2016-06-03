BINDIR = bin
SRCDIR = src
CXX = g++
CXXFLAGS = -O3 -pedantic -Wall -ansi #-static

TARGETS = correlate_dhs
EXE = $(addprefix $(BINDIR)/,$(TARGETS))

default: $(EXE)

$(BINDIR)/% : $(SRCDIR)/%.cpp
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(EXE)
