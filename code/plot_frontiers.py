import matplotlib.pyplot as plt
import sys

depth = []
size = []

filename = sys.argv[1]

file = open(filename, 'r')

runs = -1
correct_run = False
for line in file:
    if line.startswith('-------Run 5'):
        print(('-------Run '+str(runs)))
        print(line.startswith('-------Run '+str(runs)))
        print(line)
        correct_run = True
    if correct_run:
        if line.startswith('Depth:'):
            words = line.strip().split(' ')
            depth.append(int(words[1][:-1]))
            size.append(int(words[-1]))
        if line.startswith('Path length:'):
            break
    if line.startswith('Runs:'):
        runs = line.split(':')[1]
        print('looking for', '-------Run '+str(runs))
    if line.startswith('-------Run '+str(runs)):
        print('got here')
        correct_run = True

plt.figure(figsize=(9, 3))

plt.subplot(131)
plt.bar(depth, size)
plt.yscale('log')

plt.show()
