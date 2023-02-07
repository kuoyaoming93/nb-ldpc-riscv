#!/bin/bash

gcc -o example main.c -Wall -lm
rm -rf RESULTS_ETMM_q32_dc26_PRUEBA/*
./example 1977 20 3 3.1 3.2 3.3 3.4