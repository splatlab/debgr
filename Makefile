CXX = g++ -std=c++11
CC = g++ -std=c++11

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-Ofast
endif

ifdef P
	PROFILE=-pg
endif

#CXXFLAGS = -Wall -g -pg  -m64 -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function
#CXXFLAGS = -Wall -Ofast -m64 -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -DLOG_WAIT_TIME -DLOG_CLUSTER_LENGTH
CXXFLAGS = -Wall $(DEBUG) $(PROFILE) $(OPT) -m64 -I. -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function

LDFLAGS = $(DEBUG) $(PROFILE) $(OPT) -lpthread -lssl -lcrypto -lboost_system -lboost_thread -lm
#LDFLAGS = -g -pg -lpthread -lssl -lcrypto -lboost_system -lboost_thread -lm
LIBS = libs/libbz2.a libs/libz.a

TARGET_MAIN	= main
MAIN_SRC = main.cc hashutil.cc threadsafe-gqf/gqf.c

TARGET_DEBRUIJN_GRAPH = debruijn_graph
TARGET_DEBRUIJN_GRAPH_SRC = debruijn_graph.cc hashutil.cc threadsafe-gqf/gqf.c

$(TARGET_MAIN): $(MAIN_SRC)
	$(CXX) $(CXXFLAGS) $(MAIN_SRC) $(INCLUDE) $(LDFLAGS) $(LIBS) -o $@

$(TARGET_DEBRUIJN_GRAPH): $(TARGET_DEBRUIJN_GRAPH_SRC)
	$(CXX) $(CXXFLAGS) $(TARGET_DEBRUIJN_GRAPH_SRC) $(INCLUDE) $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -f $(TARGET_MAIN) $(TARGET_DEBRUIJN_GRAPH) *.o core
