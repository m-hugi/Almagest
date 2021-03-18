#!/bin/bash

gfortran -Wall -c -fdec-math almagest.f95
gfortran -Wall -o almagest main.f95 almagest.o