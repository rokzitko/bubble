CXX ?= c++
PKG_CONFIG ?= pkg-config

CXXFLAGS ?= -O2 -std=c++11
GSL_CFLAGS = $(shell $(PKG_CONFIG) --cflags gsl)
GSL_LIBS = $(shell $(PKG_CONFIG) --libs gsl)

TARGET := bubble
OBJECTS := bubble.o Faddeeva.o

.PHONY: all check clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(GSL_LIBS)

bubble.o: bubble.cc Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) -c -o $@ $<

Faddeeva.o: Faddeeva.cc Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) -c -o $@ $<

check: $(TARGET)
	@output="$$(cd regression && ./run_tests 1)"; \
	printf '%s\n' "$$output"; \
	printf '%s\n' "$$output" | grep -Fq 'OK=24 FAILED=0'
	@cd regression && ./run_optical_tests
	@cd Bethe_lattice_test && ./cond.opt.geo --check

clean:
	$(RM) $(TARGET) $(OBJECTS)
