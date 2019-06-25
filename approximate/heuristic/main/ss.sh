#!/bin/sh
cp ../../bpa16ap.blif ../verify/bpa16ap.blif
python3 main.py heu ../verify/m16c.blif ../verify/bpa16ap.blif
