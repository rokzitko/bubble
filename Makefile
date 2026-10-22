CXX ?= c++
PKG_CONFIG ?= pkg-config

CXXFLAGS ?= -O2 -std=c++11
THREAD_FLAGS ?= -pthread
GSL_CFLAGS = $(shell $(PKG_CONFIG) --cflags gsl)
GSL_LIBS = $(shell $(PKG_CONFIG) --libs gsl)

TARGET := bubble
OBJECTS := bubble.o Faddeeva.o
CLEAN_KERNEL_TEST := regression/clean_kernel_tests
FADDEEVA_TEST := regression/faddeeva_tests
LIBRARY_TEST := regression/library_tests

.PHONY: all check clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(GSL_LIBS) $(THREAD_FLAGS)

bubble.o: bubble.cc bubble.hh Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) $(THREAD_FLAGS) -c -o $@ $<

Faddeeva.o: Faddeeva.cc Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) -c -o $@ $<

bubble_test.o: bubble.cc bubble.hh Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) $(THREAD_FLAGS) -DBUBBLE_NO_MAIN -c -o $@ $<

$(CLEAN_KERNEL_TEST): regression/clean_kernel_tests.cc bubble_test.o Faddeeva.o
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(GSL_LIBS) $(THREAD_FLAGS)

$(FADDEEVA_TEST): Faddeeva.cc Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -DTEST_FADDEEVA -o $@ Faddeeva.cc $(LDLIBS) -lm

bubble_library.o: bubble.cc bubble.hh Faddeeva.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) $(THREAD_FLAGS) -DBUBBLE_LIBRARY -c -o $@ $<

$(LIBRARY_TEST): regression/library_tests.cc bubble_library.o Faddeeva.o bubble.hh
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) -I. $(LDFLAGS) -o $@ regression/library_tests.cc bubble_library.o Faddeeva.o $(LDLIBS) $(GSL_LIBS) $(THREAD_FLAGS)

check: $(TARGET) $(CLEAN_KERNEL_TEST) $(FADDEEVA_TEST) $(LIBRARY_TEST)
	@perl regression/run_library_tests
	@./$(FADDEEVA_TEST) >/dev/null
	@cd regression && perl ./check_version
	@cd regression && ./run_runner_tests
	@cd regression && ./run_tests --all
	@cd regression && ./run_optical_tests
	@cd regression && ./clean_kernel_tests
	@cd regression && ./run_input_tests
	@cd regression && perl ./run_window_tests
	@cd regression && perl ./run_dc_mesh_tests
	@cd Bethe_lattice_test && ./cond.opt.geo --check
	@cd Bethe_lattice_test && perl ./check_transport

clean:
	$(RM) $(TARGET) $(OBJECTS) bubble_test.o bubble_library.o $(CLEAN_KERNEL_TEST) $(FADDEEVA_TEST) $(LIBRARY_TEST)
