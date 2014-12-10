#!/usr/bin/python

import os;

def log_command(command, ignorefail=False, printcommand=True):
    code = os.system(command);
    if (code != 0):
        if (printcommand):
            log_print("Ran: " + command);
            log_print("Exit Code: " + str(code));
        if not ignorefail:
            raise Exception('Failed to run \'' + command + '\', exited with code ' + str(code));
    return code;

def log_print(print_str):
    os.system("echo " + str(print_str))
