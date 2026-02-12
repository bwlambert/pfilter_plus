CC=gcc
CXX=g++
RM=rm -f

CPPFLAGS= -std=c++17 -Wall -O3 -march=native -msse2 -Ieigen/ -Iboost/

LDFLAGS= -v
LDLIBS=

SRCS=pfpp.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: amain

amain: $(OBJS)
	$(CXX) $(LDFLAGS) -o pfpp $(OBJS) $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	rm pfpp.o

distclean: clean
	$(RM) *~ .depend
