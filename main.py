# main.py

import argparse
from input import (
    problem1_input,
    problem2_input,
    problem3_input,
    problem4_input,
    problem5_input
)

def main():
    parser = argparse.ArgumentParser(
        description="Run a selected problem input function with default parameters.\n" +
                    "Example: python main.py -p 4 (runs Problem 4)"
    )
    parser.add_argument(
        "-p", "--problem",
        type=int,
        choices=[1, 2, 3, 4, 5],
        required=True,
        help="Problem number to run (choose 1, 2, 3, 4, or 5)."
    )
    
    args = parser.parse_args()
    
    if args.problem == 1:
        result = problem1_input()
        print("Calculated phi for Problem 1:", result)
    elif args.problem == 2:
        result = problem2_input()
        print("Calculated phi for Problem 2:", result)
    elif args.problem == 3:
        result = problem3_input()
        print("Calculated phi for Problem 3:", result)
    elif args.problem == 4:
        result = problem4_input()
        print("Calculated phi for Problem 4:", result)
    elif args.problem == 5:
        result = problem5_input()
        print("Calculated phi for Problem 5:", result)
    else:
        print("Selected problem is not exist.")

if __name__ == "__main__":
    main()
