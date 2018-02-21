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
    # DH*
    # Nasty hack around non-CCPP-compliant module names and subroutine names
    if pgiprefix == 'radsw_cap_':
        print "WARNING! Manually change pgiprefix from 'radsw_cap_' to 'module_radsw_main_cap_'!"
        pgiprefix = 'module_radsw_main_cap_'
    elif pgiprefix == 'radlw_cap_':
        print "WARNING! Manually change pgiprefix from 'radsw_cap_' to 'module_radlw_main_cap_'!"
        pgiprefix = 'module_radlw_main_cap_'
    elif pgiprefix == 'sfc_ex_coef_cap_':
        print "WARNING! Manually change pgiprefix from 'sfc_ex_coef_cap_' to 'surface_exchange_coefficients_cap_'!"
        pgiprefix = 'surface_exchange_coefficients_cap_'
    elif pgiprefix == 'sfc_diag_cap_':
        print "WARNING! Manually change pgiprefix from 'sfc_diag_cap_' to 'surface_diagnose_cap_'!"
        pgiprefix = 'surface_diagnose_cap_'
    elif pgiprefix == 'sasasshal_cap_':
        print "WARNING! Manually change pgiprefix from 'sasasshal_cap_' to 'sasas_shal_cap_'!"
        pgiprefix = 'sasas_shal_cap_'
    elif pgiprefix == 'sasasshal_post_cap_':
        print "WARNING! Manually change pgiprefix from 'sasasshal_post_cap_' to 'sasas_shal_post_cap_'!"
        pgiprefix = 'sasas_shal_post_cap_'
    elif pgiprefix == 'gscond_cap_':
        print "WARNING! Manually change pgiprefix from 'gscond_cap_' to 'gfs_zhaocarr_gscond_cap_'!"
        pgiprefix = 'gfs_zhaocarr_gscond_cap_'
    elif pgiprefix == 'sasasdeep_cap_':
        print "WARNING! Manually change pgiprefix from 'sasasdeep_cap_' to 'sasas_deep_cap_'!"
        pgiprefix = 'sasas_deep_cap_'
    elif pgiprefix == 'precpd_cap_':
        print "WARNING! Manually change pgiprefix from 'precpd_cap_' to 'gfs_zhaocarr_precpd_cap_'!"
        pgiprefix = 'gfs_zhaocarr_precpd_cap_'
    # *DH
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