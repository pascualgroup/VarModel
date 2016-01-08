#!/usr/bin/env python

import os
import subprocess

# All paths below computed relative to location of this script using os.path.join(...)
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

def main():
    proc = subprocess.Popen([
        os.path.join(SCRIPT_DIR, 'varmodel', 'bin', 'varmodel'),
        'parameters.json'
    ])
    proc.wait()

if __name__ == '__main__':
    main()
