import itertools
import os
import subprocess

# the range of each parameter
param_range = { 'n': ['2', '10', '100', '1000', '10000', '100000', '1000000', '10000000'],
                'P': ['4', '10', '20', '30', '40'],
                't': ['0.5', '0.4', '0.3', '0.2'],
                'd': ['c', 's', 'o', 'p'] }

# the parameters of each test
test_params = { 'laplace_kernel': 'P',
                'tree': 'nd',
                'list': 'ntd',
                'laplace': 'nPtd' }

exedir_list = ['2d', '3d', '3dp']
logfile = open('nightly.log', 'w')

# loop over each directory
for exedir in exedir_list:
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
            args_list.insert(0, os.path.join(exedir, test))
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
                raise SystemExit
logfile.close()
