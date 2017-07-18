# Usage: python regression.py 3d/fmm [args]
import os
import sys
import numpy
import pickle
import subprocess
import itertools

threshold = 1 + 1e-8 + 0.01
max_error = numpy.array([5e-4, 5e-3])

def parse(lines):
    return [ [ elem.strip() for elem in line.strip().split(':') ] for line in lines ]

def accuracy_regression(args, values):
    filename = os.path.join(dirname, 'accuracy.p')
    if os.path.isfile(filename):
        print("Accuracy file exists")
        with open(filename, 'rb') as f:
            record = pickle.load(f)
    else:
        record = {}
    if (args not in record and numpy.all(values <= max_error)) or \
       (args in record and numpy.all(values <= threshold*record[args])):
        record[args] = values
        print("Accuracy regression passed")
    elif args not in record:
        print("Rel. L2 Error (p): runtime: {}, max_error: {}".format(values[0], max_error[0]))
        print("Rel. L2 Error (F): runtime: {}, max_error: {}".format(values[1], max_error[1]))
        sys.exit("Error exceeds max error")
    else:
        print("Rel. L2 Error (p): runtime: {}, record: {}".format(values[0], record[args][0]))
        print("Rel. L2 Error (F): runtime: {}, record: {}".format(values[1], record[args][1]))
        sys.exit("Accuracy regression failed")
    if len(record) > 0:
        with open(filename, 'wb') as f:
            pickle.dump(record, f)

def time_regression(args, time):
    filename = os.path.join(dirname, 'time.p')
    if os.path.isfile(filename):
        print("Time file exists")
        with open(filename, 'rb') as f:
            record = pickle.load(f)
    else:
        record = {}
    if args not in record:
        record[args] = numpy.array([time, 1])
        print("Time regression passed")
    elif time <= threshold*record[args][0]:
        record[args][0] = (record[args][0]*record[args][1] + time) / (record[args][1] + 1)
        record[args][1] += 1
        print("Time regression passed")
    else:
        print("runtime: {}, record: {}".format(time, record[args][0]))
        sys.exit("Time regression failed")
    if len(record) > 0:
        with open(filename, 'wb') as f:
            pickle.dump(record, f)

if __name__ == "__main__":
    args_list = sys.argv[1:]
    excludes = ['path', 'verbose']
    time = 0.
    repeat = 10
    dirname = os.path.dirname(sys.argv[1])
    for i in range(repeat):
        try:
            output = subprocess.check_output(args_list)
        except subprocess.CalledProcessError:
            raise SystemExit
        lines = output.strip().split("\n")
        indices = [ idx for idx, line in enumerate(lines) if line.startswith('---') ]
        # FMM Parameter section
        params = parse(lines[indices[0]+1:indices[1]])
        args = [ [key, value] for (key, value) in params if key not in excludes ]
        args = '_'.join(list(itertools.chain(*args)))
        # FMM vs. direct section
        errors = parse(lines[indices[-1]+1:])
        values = numpy.array(list(itertools.chain(*errors))[1::2], dtype=float)
        # FMM Profiling section
        for line in lines:
            if line.startswith('Total FMM'):
                time += float(line.split(':')[-1].strip())
                break
        if i == 0:
            accuracy_regression(args, values)
    time_regression(args, time/repeat)
