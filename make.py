#!/usr/bin/env python

import os
import sys
import subprocess

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

def main():
    if sys.platform.startswith('linux'):
        c_compiler = 'gcc'
        cpp_compiler = 'g++'
        flags = []
    elif sys.platform == 'darwin':
        c_compiler = 'clang'
        cpp_compiler = 'clang++'
        flags = ['-stdlib=libc++']
    else:
        sys.stderr.write('This build script only works on Linux+GCC and Mac OS X+Clang.')
        sys.exit(1)
    
    include_dirs = [
        'zppjson/src',
        'zppdb/src',
        'zppsim/src',
        'preprocessor/include',
        'libjson'
    ]
    
    src_dirs = [
        'src',
        'zppjson/src',
        'zppdb/src',
        'zppsim/src',
    ]
    
    src_files = []
    for src_dir in src_dirs:
        for root, dirs, filenames in os.walk(src_dir):
            for filename in filenames:
                if filename.endswith('.cpp'):
                    src_files.append(os.path.join(root, filename))
    
    try:
        os.makedirs(os.path.join(SCRIPT_DIR, 'bin'))
    except:
        pass
    
    run_command([
        c_compiler,
        '-O3',
        '-c', 'libjson/json.c',
        '-o', 'bin/libjson.o'
    ])
    
    run_command(
        [cpp_compiler, '-O3', '-std=c++11'] + flags +
        ['-I{}'.format(os.path.expanduser(x)) for x in include_dirs] +
        ['bin/libjson.o'] +
        src_files +
        ['-o', 'bin/varmodel', '-lsqlite3']
    )
    
    os.remove('bin/libjson.o')

def run_command(cmd_and_args):
    sys.stderr.write(' '.join(cmd_and_args))
    sys.stderr.write('\n')
    
    proc = subprocess.Popen(
        cmd_and_args,
        cwd=SCRIPT_DIR
    )
    proc.wait()

if __name__ == '__main__':
    main()
