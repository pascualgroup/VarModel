#!/usr/bin/env python

import os
import sys
import json
from collections import OrderedDict

# All paths below computed relative to location of this script using os.path.join(...)
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

def main():
    # Load existing parameters file
    params = load_json(os.path.join(SCRIPT_DIR, 'base_parameters.json'))
    
    # Information used by `runmany` to 
    runmany_info = {
        'executable' : os.path.join(SCRIPT_DIR, 'run_single.py')
    }
    
    for pMutation in [0.01, 0.02, 0.03, 0.04]:
        for immigrationRate in [0.1, 0.2, 0.3]:
            for replicate in range(1, 11):
                # Set the parameter values
                params['pMutation'] = pMutation
                params['populations'][0]['immigrationRate'] = immigrationRate
                
                # String formatting:
                # {} is a placeholder for arguments in format();
                # {:g} means general-purpose real number formatting;
                # {:d} means integer
                # {:.2g} means real to 2 decimal places;
                # {0}, {1}, {2} can be used instead of {} to specify the position of arguments;
                # and similarly {0:.2g} means use the first argument with 2-decimal-place formatting.
                
                # This example generates directories 0.01-0.05/01
                run_subdir = 'pm={:.2g}-imm={:.1g}/{:02d}'.format(pMutation, immigrationRate, replicate)
                run_dir = os.path.join(SCRIPT_DIR, 'runs', run_subdir)
                os.makedirs(run_dir)
                
                # Write out parameters
                dump_json(params, os.path.join(run_dir, 'parameters.json'))
                
                # Write out runmany info file
                dump_json(runmany_info, os.path.join(run_dir, 'runmany_info.json'))

def dump_json(obj, filename):
    with open(filename, 'w') as f:
        json.dump(obj, f, indent=4)
        f.write('\n')

def load_json(filename):
    with open(filename) as f:
        return json.load(f, object_pairs_hook=OrderedDict)

if __name__ == '__main__':
    main()
