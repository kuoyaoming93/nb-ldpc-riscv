#!/bin/bash

gcc -o example main.c -Wall -lm
rm -rf RESULTS_ETMM_q16_dc4/*
./example 1977 25 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5