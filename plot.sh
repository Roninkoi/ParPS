#!/bin/sh

# difference between serial and parallel
./plotmat.py run/serial.dat run/parallel.dat

# gather data from log files for plotting
# n, np, t
grep converged run/sinsq*.log | sed -e 's/.*sinsq//g' -e 's/n/ /' -e 's/.log:/ /g' | awk '{print $1, $2, $9}' > scaling.log
# n, gamma, it
grep gamma, run/gamma*.log | grep -v .*gamma.*n.*log.* | sed -e 's/.log:.*iterations:/ /g' -e 's/.*gamma//g' > gamma.log

# plot gamma curves
./plotgamma.py
# plot processor number scaling curves
./plotscaling.py

