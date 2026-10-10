# main.py

import argparse
from bfp.input import (
    problem1_input,
    problem2_input,
    problem3_input,
    problem4_input,
    problem5_input
)

def main():
    description=(
        "Run a selected problem input function with customizable parameters.\n\n"
        "Available problems:\n"
        "  1: Infinite Medium: ψ = Q/σₜ\n"
        "  2: Exponential Attenuation: ψ(x) = Q/σₜ + ψₗ * exp(-σₜ * x / μ)\n"
        "  3: MMS - Linear in x: ψ = a + b * x\n"
        "  4: MMS - Linear in E: ψ = a + b * E\n"
        "  5: Mixed: ψ = a + b * xE\n\n"
        "The following parameters can be provided:\n"
        "  nx (int): Number of cells in the x-direction (default: 10).\n"
        "  nE (int): Number of energy cells (default: 10).\n"
        "  N_ang (int): Number of angles for the SN method (default: 2).\n"
        "  iter_ (int): Maximum number of solver iterations (default: 1000).\n"
        "  tol (float): Solver tolerance (default: 1e-12).\n"
        "  p_level (int): Print level (1 for verbose, 0 for silent; default: 1).\n\n"
        "Example usage:\n"
        "  python main.py -p 4 --nx 15 --nE 12"
    )

    parser = argparse.ArgumentParser(
    description=description,
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-p", "--problem",
        type=int,
        choices=[1, 2, 3, 4, 5],
        required=True,
        help="Problem number to run (choose 1, 2, 3, 4, or 5)."
    )
    parser.add_argument(
        "--nx",
        type=int,
        default=10,
        help="Number of cells in the x-direction (default: 10)."
    )
    parser.add_argument(
        "--nE",
        type=int,
        default=10,
        help="Number of energy cells (default: 10)."
    )
    parser.add_argument(
        "--N_ang",
        type=int,
        default=2,
        help="Number of angles for the SN method (default: 2)."
    )
    parser.add_argument(
        "--iter_",
        type=int,
        default=1000,
        help="Maximum number of solver iterations (default: 1000)."
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-12,
        help="Solver tolerance (default: 1e-12)."
    )
    parser.add_argument(
        "--p_level",
        type=int,
        default=1,
        help="Print level (1 for verbose, 0 for silent; default: 1)."
    )
    
    args = parser.parse_args()
    problem_kwargs = dict(nx=args.nx, nE=args.nE, N_ang=args.N_ang,
                          iter_=args.iter_, tol=args.tol, p_level=args.p_level,
                          return_info=True)
    
    if args.problem == 1:
        _, info = problem1_input(**problem_kwargs)
    elif args.problem == 2:
        _, info = problem2_input(**problem_kwargs)
    elif args.problem == 3:
        _, info = problem3_input(**problem_kwargs)
    elif args.problem == 4:
        _, info = problem4_input(**problem_kwargs)
    elif args.problem == 5:
        _, info = problem5_input(**problem_kwargs)
    else:
        print("Selected problem is not exist.")
        return

    failed = [entry['mu'] for entry in info if not entry['converged']]
    if failed:
        print(f"Problem {args.problem}: GMRES did not converge for {len(failed)} "
              "angular direction(s). Saved results are unconverged.")
    else:
        print(f"Problem {args.problem}: GMRES converged for all angular directions. Results saved.")

if __name__ == "__main__":
    main()
