#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: plot.py file1.dat [file2.dat ...]")
    sys.exit(1)

for path in sys.argv[1:]:
    data = np.loadtxt(path)
    plt.plot(data[:, 0], data[:, 1], label=path)

plt.legend()
plt.tight_layout()
plt.show()
