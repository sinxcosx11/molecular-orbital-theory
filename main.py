# mo_diagram_generator/main.py

# (Imports remain the same)
import sys
import os
import traceback
from typing import Optional # Added for type hinting

try:
    from .cli import setup_arg_parser
    from .utils import parse_molecule_input
    # Import new function get_neutral_valence_electrons
    from .atoms import get_atom_data, get_total_valence_electrons, get_neutral_valence_electrons
    from .orbitals import generate_molecular_orbitals
    from .electron_filler import fill_electrons
    from . import properties
    from .plotter import plot_mo_diagram
except ImportError as e:
     # (Error handling remains the same)
     if "__main__" == __name__ and not __package__:
          script_path = os.path.abspath(__file__)
          parent_dir = os.path.dirname(os.path.dirname(script_path))
          module_name = os.path.basename(os.path.dirname(script_path))
          print(f"Import Error: {e}.\n\nIt seems you are running '{os.path.basename(script_path)}' directly.\n"
                f"Navigate one directory up (to '{parent_dir}') and use:\n"
                f"  python -m {module_name}\n\n"
                f"The script will prompt you for the molecule name.\n", file=sys.stderr)
     else:
          print(f"Import Error: {e}. Check package structure and dependencies.", file=sys.stderr)
     sys.exit(1)


def run_analysis(molecule_input_str: str, output_filename_arg: Optional[str]):
    """Performs the MO analysis for a given molecule string."""
    try:
        # 1. Parse Input (remains the same)
        atom1_sym, atom2_sym, charge = parse_molecule_input(molecule_input_str)
        is_homonuclear = (atom1_sym == atom2_sym)
        charge_str = ""
        if charge == 1: charge_str = "+"
        elif charge == -1: charge_str = "-"
        elif charge > 1: charge_str = f"{charge}+"
        elif charge < -1: charge_str = f"{-charge}-"
        molecule_label_for_title = f"{atom1_sym}{atom2_sym}{charge_str}"
        original_input_label = molecule_input_str

        print(f"\n--- Analyzing Species: {molecule_label_for_title} (from input '{original_input_label}') ---")
        print(f"   Type: {'Homonuclear' if is_homonuclear else 'Heteronuclear'} Diatomic")
        print(f"   Charge: {charge:+}")

        # 2. Get Atomic Data (remains the same)
        print("\nFetching atomic data...")
        atom_A_data = get_atom_data(atom1_sym)
        atom_B_data = get_atom_data(atom2_sym)
        print(f"Atom A ({atom1_sym}): Z={atom_A_data['atomic_number']}, Period={atom_A_data['period']}, Group={atom_A_data['group']}, Valence Shells: {atom_A_data['valence_orbital_types']}")
        print(f"Atom B ({atom2_sym}): Z={atom_B_data['atomic_number']}, Period={atom_B_data['period']}, Group={atom_B_data['group']}, Valence Shells: {atom_B_data['valence_orbital_types']}")

        # 3. Calculate Valence Electrons (Get counts for plotting)
        # Get neutral counts FIRST for AO plotting
        ve_A_neutral = get_neutral_valence_electrons(atom_A_data)
        ve_B_neutral = get_neutral_valence_electrons(atom_B_data)
        # Then calculate total VE for the species (including charge) for MO filling
        total_valence_electrons = get_total_valence_electrons(atom_A_data, atom_B_data, charge)
        print(f"\nNeutral Atom Valence Electrons: {atom1_sym}={ve_A_neutral}, {atom2_sym}={ve_B_neutral}")
        print(f"Total Valence Electrons to fill in MOs (Charge {charge:+}): {total_valence_electrons}")


        # 4. Generate Atomic and Molecular Orbitals (remains the same)
        print("\nGenerating Atomic and Molecular Orbitals...")
        all_original_aos, mos = generate_molecular_orbitals(atom_A_data, atom_B_data)

        if not mos:
             print("\nError: No molecular orbitals were generated.", file=sys.stderr)
             return False

        # 5. Fill Electrons (remains the same)
        print(f"\nFilling {total_valence_electrons} electrons into MOs...")
        electron_config_str = fill_electrons(mos, total_valence_electrons)
        print(f"\nElectron Configuration:\n  {electron_config_str}")

        # 6. Calculate Properties (remains the same)
        print(f"\nCalculated Properties:")
        bond_order = properties.calculate_bond_order(mos)
        magnetic_behavior = properties.determine_magnetic_behavior(mos)
        homo, lumo = properties.find_homo_lumo(mos)
        homo_lumo_gap = properties.calculate_homo_lumo_gap(homo, lumo)

        print(f"  - Bond Order: {bond_order:.2f}")
        print(f"  - Magnetic Behavior: {magnetic_behavior}")
        if homo: print(f"  - HOMO: {homo.symbol} (Energy: {homo.energy:.3f} eV)")
        else: print("  - HOMO: Not found")
        if lumo: print(f"  - LUMO: {lumo.symbol} (Energy: {lumo.energy:.3f} eV)")
        else: print("  - LUMO: Not found")
        if homo_lumo_gap is not None: print(f"  - HOMO-LUMO Gap: {homo_lumo_gap:.3f} eV")
        else: print("  - HOMO-LUMO Gap: Not applicable")

        # 7. Plot the MO Diagram (Pass neutral VE counts)
        output_filename = output_filename_arg
        if not output_filename:
            # (Filename generation logic remains the same)
            safe_mol_name = "".join(c if c.isalnum() else '_' for c in molecule_label_for_title)
            safe_mol_name = safe_mol_name.replace('+','_plus').replace('-','_minus')
            output_filename = f"{safe_mol_name}_MO_diagram.svg"

        # (Filename validation logic remains the same)
        supported_formats = ['.png', '.svg', '.pdf', '.jpg', '.jpeg']
        try:
             file_base, file_ext = os.path.splitext(output_filename)
             if not file_ext or file_ext.lower() not in supported_formats:
                  original_name = output_filename
                  output_filename = f"{file_base}.svg"
                  print(f"\nWarning: Output filename '{original_name}' lacked supported extension. Using '{output_filename}'.")
        except Exception:
             print(f"Warning: Could not parse output filename '{output_filename}'. Using default.", file=sys.stderr)
             safe_mol_name = "".join(c if c.isalnum() else '_' for c in molecule_label_for_title)
             safe_mol_name = safe_mol_name.replace('+','_plus').replace('-','_minus')
             output_filename = f"{safe_mol_name}_MO_diagram.svg"

        # Filter AOs for plotting (remains the same)
        aos_A_for_plot = [ao for ao in all_original_aos if ao.atom_label == 'A']
        aos_B_for_plot = [ao for ao in all_original_aos if ao.atom_label == 'B']

        print(f"\nGenerating plot: '{output_filename}'...")
        # ---> PASS ve_A_neutral and ve_B_neutral to plotter <---
        plot_mo_diagram(
            aos_A=aos_A_for_plot,
            aos_B=aos_B_for_plot,
            ve_A=ve_A_neutral, # Pass neutral valence electron count for Atom A
            ve_B=ve_B_neutral, # Pass neutral valence electron count for Atom B
            mos=mos,
            atom_A_label=atom1_sym,
            atom_B_label=atom2_sym,
            title=molecule_label_for_title,
            bond_order=bond_order,
            magnetic_behavior=magnetic_behavior,
            output_filename=output_filename
        )

        print("\n--- Analysis Complete ---")
        return True

    # (Exception handling remains the same)
    except ValueError as e:
        print(f"\nInput or Calculation Error: {e}", file=sys.stderr)
    except FileNotFoundError as e:
         print(f"\nFile Error: {e}. Check required data files.", file=sys.stderr)
    except Exception as e:
        print(f"\nAn unexpected error occurred during analysis: {e}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        print("Please check input and data.", file=sys.stderr)

    return False


def run_interactive():
    """Main interactive execution loop."""
    # (Interactive loop logic remains the same)
    parser = setup_arg_parser()
    args = parser.parse_args()

    while True:
        try:
            molecule_input = input("\nEnter the diatomic molecule or ion (e.g., O2, N2-, CO) or type 'quit' to exit: ")
            molecule_input = molecule_input.strip()
            if not molecule_input: continue
            if molecule_input.lower() == 'quit': break

            success = run_analysis(molecule_input, args.output)
            if not success:
                 print(f"\nAnalysis failed for '{molecule_input}'. Please check the input format and try again.")

        except KeyboardInterrupt:
            print("\nExiting due to user interrupt.")
            break
        except EOFError:
             print("\nExiting.")
             break
        except Exception as e:
            print(f"\nAn unexpected error occurred in the main loop: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)


# --- Entry point ---
if __name__ == "__main__":
    run_interactive()