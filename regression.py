# Usage: python regression.py [dir]/fmm [args]
import os
import sys
import numpy
import pickle
import subprocess
import itertools

def parse(lines):
    return [ [ elem.strip() for elem in line.strip().split(':') ] for line in lines ]

def print_helper(label, runtime, benchmark_type, benchmark):
    '''label (str): potential / force
       benchmark_type (str): max_error or record'''
    print("{}: runtime: {}, {}: {}".format(label, runtime, benchmark_type, benchmark))

def accuracy_regression(kernel, args, labels, values, max_error):
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
        for i in xrange(len(labels)):
            print_helper(labels[i], values[i], 'max_error', max_error[i])
        print("Accuracy regression passed")
    elif args not in record:
        for i in xrange(len(labels)):
            print_helper(labels[i], values[i], 'max_error', max_error[i])
        sys.exit("Error exceeds max error")
    else:
        for i in xrange(len(labels)):
            print_helper(labels[i], values[i], 'record', record[args][i])
        sys.exit("Accuracy regression failed")
    if len(record) > 0:
        with open(filename, 'wb') as f:
            pickle.dump(record, f)

if __name__ == "__main__":
    args_list = sys.argv[1:]
    excludes = ['path', 'verbose']
    dirname = os.path.dirname(sys.argv[0])
    try:
        output = subprocess.check_output(args_list)
    except subprocess.CalledProcessError:
        raise SystemExit
    lines = output.strip().split("\n")
    indices = [ idx for idx, line in enumerate(lines) if line.startswith('---') ]
    # FMM Parameter section
    params = parse(lines[indices[0]+1:indices[1]])
    max_error = numpy.array([5e-4, 5e-5])
    if map(int, [ value for (key, value) in params if key in 'numBodies']) < [10]:
        max_error = numpy.array([1e-3, 1e-2])
    args = [ [key, value] for (key, value) in params if key not in excludes ]
    if 'mpirun' in sys.argv:
        idx = sys.argv.index('mpirun')
        nproc = ['nproc', sys.argv[idx+2]]
        args.insert(0, nproc)
        del sys.argv[idx:idx+3]
    kernel = os.path.split(sys.argv[1])[-1]
    args = kernel + '_' + '_'.join(list(itertools.chain(*args)))
    # FMM vs. direct section
    errors = list(itertools.chain(*parse(lines[indices[-1]+1:])))
    labels = errors[::2]
    values = numpy.array(errors[1::2], dtype=float)
    accuracy_regression(kernel, args, labels, values, max_error)
