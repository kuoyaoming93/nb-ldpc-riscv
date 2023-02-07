#!/bin/bash

gcc -o example main.c -Wall -lm
rm -rf RESULTS_ETMM_q32_dc26_NUEVO/*
./example 1977 25 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5