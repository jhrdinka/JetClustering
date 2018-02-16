PATH := $(PATH):$(PWD)/fastjet/bin
SHELL := env PATH=$(PATH) /bin/bash

CONTRIB_INC=-I${PWD}/fastjet/include/fastjet/contrib/
CONTRIB_LIB=-I${PWD}/fastjet/lib -lRecursiveTools -lNsubjettiness -lSoftKiller

CLASSES_INC=-I${PWD}/classes

CC=gcc
CXX=g++
RM=rm -f

CPPFLAGS=-g $(shell root-config --cflags) $(shell fastjet-config --cxxflags) $(CONTRIB_INC) $(CLASSES_INC) -std=c++1y -Wall
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=-g $(shell root-config --libs) $(shell fastjet-config --libs) $(CONTRIB_LIB)

SRCS=$(wildcard src/*.cc) $(wildcard classes/*.cc) 
OBJS=$(subst .cc,.o,$(SRCS))

all: analyze analyzeCluster

analyze: $(OBJS)
	$(CXX) $(LDFLAGS) -o analyze $(OBJS) $(LDLIBS) 

analyzeCluster: $(OBJS)
	$(CXX) $(LDFLAGS) -o analyzeCluster $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend
