# Usage: python regression.py [dir]/fmm [args]
import os
import sys
import numpy
import pickle
import subprocess
import itertools

def parse(lines):
    return [ [ elem.strip() for elem in line.strip().split(':') ] for line in lines ]

def accuracy_regression(args, values, max_error):
    filename = os.path.join(dirname, 'accuracy.p')
    if os.path.isfile(filename):
        print("Accuracy file exists")
        with open(filename, 'rb') as f:
            record = pickle.load(f)
    else:
        record = {}
    if (args not in record and numpy.all(values <= max_error)) or \
       (args in record and numpy.all(values <= 1.000001*record[args]+1e-8)):
        record[args] = values
        print("Rel. L2 Error (p): runtime: {}, max_error: {}".format(values[0], max_error[0]))
        print("Rel. L2 Error (F): runtime: {}, max_error: {}".format(values[1], max_error[1]))
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

if __name__ == "__main__":
    args_list = sys.argv[1:]
    excludes = ['path', 'verbose']
    dirname = os.path.dirname(sys.argv[1])
    try:
        output = subprocess.check_output(args_list)
    except subprocess.CalledProcessError:
        raise SystemExit
    lines = output.strip().split("\n")
    indices = [ idx for idx, line in enumerate(lines) if line.startswith('---') ]
    # FMM Parameter section
    params = parse(lines[indices[0]+1:indices[1]])
    max_error = numpy.array([6e-5, 2e-5])
    if map(int, [ value for (key, value) in params if key in 'numBodies']) < [10]:
        max_error = numpy.array([1e-3, 1e-2])
    args = [ [key, value] for (key, value) in params if key not in excludes ]
    args = '_'.join(list(itertools.chain(*args)))
    # FMM vs. direct section
    errors = parse(lines[indices[-1]+1:])
    values = numpy.array(list(itertools.chain(*errors))[1::2], dtype=float)
    accuracy_regression(args, values, max_error)
