#!/bin/bash

gcc -o example main.c -Wall -lm
rm -rf RESULTS_ETMM_q16_dc4_float/*
./example 1977 25 5