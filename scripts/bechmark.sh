#!/bin/bash
hostname
strings Benchmark | grep " -m"

for lineSize in 1024 512 256 128 64; do
  for offsetSize in 32 16; do
    ./Benchmark --numObjects 100M --numQueries 100M --lineSize $lineSize --offsetSize $offsetSize
  done
done
