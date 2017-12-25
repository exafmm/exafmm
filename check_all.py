import itertools
import os
import subprocess
from collections import OrderedDict

# the range of each parameter
param_range = { 'n': ['2', '10', '100', '1000'],
                'P': ['20'],
                't': ['0.5', '0.4'],
                'd': ['c', 's', 'p', 'l'] }

# the parameters of each test
test_params = OrderedDict([('laplace_kernel', 'P'),
                           ('laplace_ki_kernel', 'P'),
                           ('helmholtz_kernel', 'P'),
                           ('stokes_kernel', 'P'),
                           ('build_tree', 'nd'),
                           ('traverse', 'ntd'),
                           ('laplace', 'nPtd'),
                           ('laplace_ki', 'nPtd'),
                           ('helmholtz', 'nPtd'),
                           ('stokes', 'nPtd')])

testdir = "tests"
logfile = open('check_all.log', 'w')

# loop over each test
for test, params in test_params.iteritems():
    # concatenate ranges for itertools
    ranges_list = [ param_range[param] for param in params]
    # nested loop of each parameter's range
    for values in itertools.product(*ranges_list):
        # add dash before each parameter character (ex. 'n'->'-n')
        dash_params = [ '-'+param for param in params]
        # interleave parameter character list and their value list
        args_list = list(itertools.chain(*zip(dash_params, values)))
        # insert executable at the beginning of args list
        args_list.insert(0, os.path.join(testdir, test))
        args_list.insert(0, "4")
        args_list.insert(0, "-np")
        args_list.insert(0, "mpirun")
        # for fmm executable
        if test is "fmm" or test is "laplace":
            args_list.insert(0, "regression.py")
            args_list.insert(0, "python")
        # print args string and write to log file
        args_str = ' '.join(args_list)
        print(args_str)
        logfile.write(args_str + 2*'\n')
        # execute test
        try:
            result = subprocess.check_output(args_list)
            logfile.write(result)
            logfile.write(33*'-' + '\n')
        except subprocess.CalledProcessError:
            print("Regression Failed @", args_str)
            print("...trying again...")
            subprocess.call(args_list)
            raise SystemExit
logfile.close()
