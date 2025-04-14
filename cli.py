# mo_diagram_generator/cli.py
import argparse
import os

def setup_arg_parser():
    """Sets up the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate Molecular Orbital diagrams for diatomic species interactively.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Show defaults in help
        epilog="Example: Run 'python -m mo' and follow the prompts." # Updated epilog
    )
    # REMOVED: The positional 'molecule' argument is no longer needed here.
    # parser.add_argument(
    #     "molecule",
    #     type=str,
    #     help="The diatomic molecule or ion to analyze (e.g., 'O2', 'N2-', 'CO', 'LiH', 'CN+'). "
    #          "For single atom symbols like 'N-', it assumes 'N2-'."
    # )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None, # Default filename generated dynamically in main.py
        help="Optional: Output filename for the diagram (e.g., 'my_diagram.png'). "
             "Supported formats: png, svg, pdf, jpg. If not specified or invalid, "
             "defaults to '<molecule_name>_MO_diagram.svg'."
    )
    # Potential future arguments can remain commented out or be added here if needed.
    # parser.add_argument(
    #     "--beta",
    #     type=float,
    #     # default=orbitals.DEFAULT_BETA, # Needs import or passing value
    #     help="Interaction parameter (beta) in eV, controlling MO splitting."
    # )
    # parser.add_argument(
    #     "--energy-threshold",
    #     type=float,
    #     # default=orbitals.ENERGY_DIFF_THRESHOLD, # Needs import or passing value
    #     help="Energy difference threshold (eV) for AO interaction."
    # )

    return parser