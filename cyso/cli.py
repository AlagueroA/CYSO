# src/cli.py
import argparse
import shutil
import sys
from pathlib import Path
import importlib.resources as resources



def main():

    print('')
    print('')
    print('')
    print('')
    print(' ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓███████▓▒░░▒▓██████▓▒░')
    print('░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░')
    print('░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░')
    print('░▒▓█▓▒░       ░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░')
    print('░▒▓█▓▒░         ░▒▓█▓▒░          ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░')
    print('░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░          ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░')
    print(' ░▒▓██████▓▒░   ░▒▓█▓▒░   ░▒▓███████▓▒░ ░▒▓██████▓▒░')
    #print('Create Your Synthetic Observation')
    #print('A Python package by Antoine Alaguero')
    print('')
    print('')
    print('')
    print('')
    
    parser = argparse.ArgumentParser(prog="cyso")
    parser.add_argument("--make_setup", action="store_true",
                        help="Copy CYSO_input.py to current directory")
    parser.add_argument("--cook", action="store_true",
                        help="Run main.py using the copied CYSO_input.py")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        return
        
    if args.make_setup:
        copy_user_input()

    if args.cook:
        from cyso.main import run_main
        run_main()


def copy_user_input():
    """Copy CYSO_input.py from the package to the current directory."""
    target = Path("CYSO_input.py")

    if target.exists():
        print("CYSO_input.py already exists in this directory. Please rm it if you want to create a new one.")
        return

    # Access the packaged file using importlib.resources
    with resources.path("cyso", "CYSO_input.py") as src:
        shutil.copy(src, target)

    print("CYSO_input.py has been created in the current directory. Modify it to your convenience and run cyso --cook.")
