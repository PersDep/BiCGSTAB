#!/usr/bin/python

from collections import defaultdict
import matplotlib.pyplot as plotter

data = open("big_result.out", "r").readlines()

table = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

threads_amounts = [1, 2, 4, 8]
figures = ['_', '+', 'x', '|']
names = ['Global time', 'Linear combination time', 'Scalar produc time', 'Mutiplication time']

index = -1
size = -1
threads = 0
for line in data:
    if 'Iterations' in line:
        index = -1
    if index != -1:
        cols = line.split(':')
        if len(cols) != 2:
            cols = line.split(' ')
            size = int(cols[0]) * int(cols[1]) * int(cols[2])
            continue
        if cols[0] == "Threads":
            threads = int(cols[1])
            continue
        table[size][index][threads] = float(cols[1])
        index += 1
    if 'Start' in line:
        index = 0

for size, types in sorted(table.items()):
    for type, threads in types.items():
        by_threads = []
        for threads_amount, time in sorted(threads.items()):
            print(str(size) + ': ' + str(type) + ' ' + str(threads_amount) + ' ' + str(time))
            by_threads.append(time)
        plotter.xlabel('Threads amount')
        plotter.ylabel('Time for size ' + str(size) + ', sec')
        plotter.grid('True')
        plotter.plot(threads_amounts, by_threads, 'c', linewidth = '0.5')
        plotter.plot(threads_amounts, by_threads, 'b' + figures[type], label=names[type], markersize='15')
    plotter.legend()
    plotter.show()

