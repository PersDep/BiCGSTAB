#!/usr/bin/python

from collections import defaultdict
import matplotlib.pyplot as plotter

time_table = defaultdict(lambda: defaultdict(list))
acceleration_table = defaultdict(lambda: defaultdict(list))

sizes = [100000, 1000000, 10000000]
nodes_amount = [1, 2, 4, 8, 16]
acceleration_names = ['Global acceleration', 'Linear combination acceleration', 'Scalar product acceleration', 'Mutiplication acceleration']
figures = ['_', '+', 'x', '|']
colors = ['g', 'r', 'c']

output = False

for size in sizes:
    for i in nodes_amount:
        data = open("result_mpi_" + str(i) + '_' + str(size) + ".out", "r").readlines()
        time_type = 0
        for line in data:
            if 'Time' in line or 'time' in line:
                columns = line.split(':')
                time = float(columns[1])
                time_table[size][time_type].append(time)
                time_type += 1

for size, types in sorted(time_table.items()):
    for type, times in sorted(types.items()):
        for time in times:
            acceleration_table[size][type].append(times[0] / time)

if output:
    for size, types in sorted(time_table.items()):
        for type, times in sorted(types.items()):
            print(size, type)
            print(times)
            print()
    for size, types in sorted(acceleration_table.items()):
        for type, acceleration in sorted(types.items()):
            print(size, type)
            print(acceleration)
            print()

acceleration_ticks = []
color_index = -1
for size, types in sorted(acceleration_table.items()):
    color_index += 1
    for type, acceleration in sorted(types.items()):
        plotter.xlabel('Nodes amount')
        plotter.ylabel('Acceleration')
        plotter.xscale('log')
        plotter.grid(True)
        acceleration_ticks += acceleration
        plotter.plot(nodes_amount, acceleration, colors[color_index], linewidth='1')
        if type == 0:
            plotter.plot(nodes_amount, acceleration, colors[color_index], label="Size: " + str(size), linewidth='1')
        plotter.plot(nodes_amount, acceleration, 'b' + figures[type], markersize='15')
        if size == sizes[2]:
            plotter.plot(nodes_amount, acceleration, 'b' + figures[type], label=acceleration_names[type], markersize='15')
            plotter.text(nodes_amount[0] - 0.05, acceleration[0] - 0.5, 'time for size 10^7 on 1 node is extrapolated', color=colors[color_index])

ticks_to_remove = []
acceleration_ticks.sort()
for i, tick in enumerate(acceleration_ticks):
    if i > 0 and tick - acceleration_ticks[i - 1] < 0.2:
        ticks_to_remove.append(tick)
for tick in ticks_to_remove:
    acceleration_ticks.remove(tick)

plotter.xticks(nodes_amount, nodes_amount)
plotter.yticks(acceleration_ticks, acceleration_ticks)
plotter.legend()
plotter.show()