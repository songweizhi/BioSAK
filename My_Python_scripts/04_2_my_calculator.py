import sys
import argparse
from my_calculator_modules.get_sum import get_sum
from my_calculator_modules.get_product import get_product
from my_calculator_modules.get_square import get_square
from my_calculator_modules.get_cube import get_cube


def print_main_help():

    help_message = '''    
    .....::::: my_calculator :::::.....
               
    add       -> calculate sum 
    multiply  -> calculate product
    square    -> calculate square
    cube      -> calculate cube
    
    Usage: my_calculator <command> -h for command specific help
    '''

    print(help_message)


# initialize the options parser
parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help="--", dest='subparser_name')


# arguments for add
add_parser = subparsers.add_parser('add',
                                   description='calculate sum of specified two numbers',
                                   epilog='Example: my_calculator add -a1 5 -a2 7')

add_parser.add_argument('-a1', required=True, type=int, help='input number 1')
add_parser.add_argument('-a2', required=True, type=int, help='input number 2')


# arguments for multiply
multiply_parser = subparsers.add_parser('multiply',
                                        description='calculate multiply',
                                        epilog='Example: my_calculator square -m1 5 -m2 6')

multiply_parser.add_argument('-m1', required=True, type=int, help='input number 1')
multiply_parser.add_argument('-m2', required=True, type=int, help='input number 2')


# arguments for square
square_parser = subparsers.add_parser('square',
                                      description='calculate square',
                                      epilog='Example: my_calculator square -s 6')

square_parser.add_argument('-s', required=True, type=int, help='input number')


# arguments for cube
cube_parser = subparsers.add_parser('cube',
                                    description='calculate cube',
                                    epilog='Example: my_calculator cube -c 5')

cube_parser.add_argument('-c', required=True, type=int, help='input number')


# get and check options
args = None
if (len(sys.argv) == 1) or (sys.argv[1] == '-h') or (sys.argv == '--help'):
    print_main_help()
    sys.exit(0)
else:
    args = parser.parse_args()
args_vars = vars(args)


if args_vars['subparser_name'] == 'add':
    get_sum(args_vars)

if args_vars['subparser_name'] == 'multiply':
    get_product(args_vars)

if args_vars['subparser_name'] == 'square':
    get_square(args_vars)

if args_vars['subparser_name'] == 'cube':
    get_cube(args_vars)

