#!/bin/bash

gfortran -Wall -c -fdec-math almagest.f03
gfortran -Wall -o almagest main.f03 almagest.o
