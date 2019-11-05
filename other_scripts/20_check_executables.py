import os
import subprocess


def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:

        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


check_executables(['mafft', 'hmmsearch'])
