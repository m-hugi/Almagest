#!/bin/bash

gfortran -c degfun.f95
gfortran -o almagest main.f95 degfun.o
