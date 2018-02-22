#!/usr/bin/env python

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='Fix cap objects produced by PGI compiler')
parser.add_argument("cap")

def parse_args():
    args = parser.parse_args()
    cap = args.cap
    return cap

def execute(cmd, debug = True, abort = True):
    """Runs a local command in a shell. Waits for completion and
    returns status, stdout and stderr. If abort = True, abort in
    case an error occurs during the execution of the command."""
    
    if debug:
        print 'Executing "{0}"'.format(cmd)
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE, shell = True)
    (stdout, stderr) = p.communicate()
    status = p.returncode
    if debug:
        message = 'Execution of "{0}" returned with exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.rstrip('\n'))
        print message
    if not status == 0:
        message = 'Execution of command {0} failed, exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.rstrip('\n'))
        if abort:
            raise Exception(message)
        else:
            print message
    return (status, stdout.rstrip('\n'), stderr.rstrip('\n'))

def correct_cap_object_names(fixcmd, cap):
    (cappath, capname) = os.path.split(cap)
    pgiprefix = capname.rstrip('.o').lower() + '_'
    # Get list of all symbols in cap object
    nmcmd = 'nm {0}'.format(cap)
    (status, stdout, stderr) = execute(nmcmd)
    del nmcmd
    # Parse all symbols and generate objcopy command
    found = False
    for line in stdout.split('\n'):
        try:
            (address, symboltype, objectname) = line.split()
        except ValueError:
            continue
        if not symboltype == 'T':
            continue
        if objectname.startswith(pgiprefix):
            newname = objectname[len(pgiprefix):]
        else:
            continue
        if newname.endswith('_cap'):
            fixcmd += '--redefine-sym {0}={1} '.format(objectname, newname)
            found = True
    if not found:
        raise Exception('Unable to rename CCPP scheme caps in cap "{0}"'.format(cap))
    return fixcmd

def correct_object_names(fixcmd, cap):
    tmp = cap + '.tmp'
    fixcmd += '{0} {1}'.format(cap, tmp)
    execute(fixcmd)
    mvcmd = 'mv -v {0} {1}'.format(tmp, cap)
    execute(mvcmd)

def main():
    cap = parse_args()
    fixcmd = 'objcopy '
    fixcmd = correct_cap_object_names(fixcmd, cap)
    if not fixcmd == 'objcopy ':
        correct_object_names(fixcmd, cap)

if __name__ == '__main__':
    main()