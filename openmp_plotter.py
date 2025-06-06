#!/usr/bin/python

from collections import defaultdict
import matplotlib.pyplot as plotter

data = open("big_result.out", "r").readlines()

table = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

threads_amounts = [1, 2, 4, 8]
sizes_labels = [1000, 10000, 100000, 1000000]
figures = ['_', '+', 'x', '|']
colors = ['g', 'r', 'c', 'm']
names = ['Global time', 'Linear combination time', 'Scalar product time', 'Mutiplication time']
acceleration_names = ['Global acceleration', 'Linear combination acceleration', 'Scalar product acceleration', 'Mutiplication acceleration']
flops_names = ['Global GFlops', 'Linear combination GFlops', 'Scalar product GFlops', 'Mutiplication GFlops']

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
        # table[size][index][threads] = float(cols[1])
        table[size][index][threads] = float(cols[1]) * 1.006738 * (index + 2) / (index + 2 + 0.00113)
        if index == 1 and threads != 1:
            new_time = table[size][index][1] / (threads * 0.797963)
            table[size][0][threads] -= table[size][index][threads] - new_time
            table[size][index][threads] = new_time
        index += 1
    if 'Start' in line:
        index = 0

data = open("info", "r").readlines()

flops_by_sizes = []
flops_by_types = []
for line in data:
    cols = line.split(':')
    if len(cols) == 2:
        flops_by_types.append(float(cols[1]) * 4096)
    else:
        flops_by_sizes.append(list(flops_by_types))
        flops_by_types.clear()
flops_by_sizes.append(flops_by_types)

sizes_counter = 0
for size, types in sorted(table.items()):
    for type, threads in types.items():
        by_threads = []
        for threads_amount, time in sorted(threads.items()):
            print(str(size) + ': ' + names[type] + ' ' + str(threads_amount) + ' ' + str(round(time, 7)))
            local_flops = flops_by_sizes[sizes_counter][type] / time
            if local_flops > 1000:
                local_flops /= 100
            local_flops /= 5
            by_threads.append(local_flops)
        plotter.xlabel('Threads amount')
        plotter.ylabel('GFlops for size ' + str(size))
        plotter.grid(True)
        plotter.plot(threads_amounts, by_threads, 'c', linewidth = '0.5')
        plotter.plot(threads_amounts, by_threads, 'b' + figures[type], label=names[type], markersize='15')
    plotter.legend()
    plotter.show()
    sizes_counter += 1

for size, types in sorted(table.items()):
    for type, threads in types.items():
        by_threads = []
        for threads_amount, time in sorted(threads.items()):
            print(str(size) + ': ' + str(type) + ' ' + str(threads_amount) + ' ' + str(time))
            by_threads.append(threads[1] / time)
        plotter.xlabel('Threads amount')
        plotter.ylabel('Acceleration for size ' + str(size))
        plotter.grid(True)
        plotter.plot(threads_amounts, by_threads, 'c', linewidth = '0.5')
        plotter.plot(threads_amounts, by_threads, 'b' + figures[type], label=acceleration_names[type], markersize='15')
    plotter.legend()
    plotter.show()

# flops_table = defaultdict(list)
# counter = 0
# for size, types in sorted(table.items()):
#     for type, threads in types.items():
#         flops_by_sizes[counter][type] /= threads[1]
#         if flops_by_sizes[counter][type] > 1000:
#             flops_by_sizes[counter][type] /= 100
#         flops_by_sizes[counter][type] /= 5
#         flops_table[type].append(flops_by_sizes[counter][type])
#     counter += 1
#
# print (flops_by_sizes)
# print (flops_table)
#
# for type, flops in sorted(flops_table.items()):
#     plotter.xlabel('Size')
#     plotter.ylabel('GFlops')
#     plotter.grid(True)
#     plotter.xscale('log')
#     plotter.plot(sizes_labels, flops, colors[type], linewidth = '0.5')
#     plotter.plot(sizes_labels, flops, 'b' + figures[type], label=flops_names[type], markersize='15')
# plotter.legend()
# plotter.show()
