# mo_diagram_generator/properties.py
from .orbitals import MolecularOrbital
from typing import List, Tuple, Optional
import sys # Added for warnings/errors if needed

def calculate_bond_order(mos: List[MolecularOrbital]) -> float:
    """
    Calculates the bond order from the list of filled MOs.
    Non-bonding electrons do not contribute.
    """
    bonding_electrons = 0
    antibonding_electrons = 0

    for mo in mos:
        if mo.is_bonding is True: # Explicitly bonding
            bonding_electrons += mo.electron_count
        elif mo.is_bonding is False: # Explicitly antibonding
            antibonding_electrons += mo.electron_count
        # Non-bonding (is_bonding is None) are ignored

    bond_order = (bonding_electrons - antibonding_electrons) / 2.0

    # Bond order cannot be negative
    return max(0.0, bond_order)


def determine_magnetic_behavior(mos: List[MolecularOrbital]) -> str:
    """
    Determines if the molecule is paramagnetic (unpaired electrons) or
    diamagnetic (all electrons paired) based on the filled MOs.
    """
    total_unpaired_electrons = 0
    for mo in mos:
        # An MO contributes to paramagnetism if it is singly occupied (has exactly 1 electron).
        if mo.electron_count == 1:
            total_unpaired_electrons += 1
        # Assumes capacity 2 orbitals. If capacity > 2 were implemented,
        # this check would need refinement based on Hund's rule application.

    if total_unpaired_electrons > 0:
        # Use Unicode minus sign for electron charge symbol
        plural = "s" if total_unpaired_electrons > 1 else ""
        return f"Paramagnetic ({total_unpaired_electrons} unpaired e⁻)"
    else:
        return "Diamagnetic"

# --- NEW FUNCTIONS FOR HOMO/LUMO ---

def find_homo_lumo(filled_mos: List[MolecularOrbital]) -> Tuple[Optional[MolecularOrbital], Optional[MolecularOrbital]]:
    """
    Finds the Highest Occupied Molecular Orbital (HOMO) and
    Lowest Unoccupied Molecular Orbital (LUMO) from a list of MOs
    ASSUMED TO BE SORTED BY ENERGY.

    Args:
        filled_mos: List of MolecularOrbital objects, sorted by energy,
                    with electron counts populated.

    Returns:
        A tuple containing (homo, lumo). Either can be None if not found.
    """
    homo: Optional[MolecularOrbital] = None
    lumo: Optional[MolecularOrbital] = None

    # Iterate through MOs sorted by energy
    for mo in filled_mos:
        # Find the highest energy orbital that is occupied
        if mo.electron_count > 0:
            # If homo is None or current MO energy is higher than current homo
            # (Handles degenerate HOMOs - takes the last one in the sorted list)
            if homo is None or mo.energy >= homo.energy:
                 homo = mo

        # Find the lowest energy orbital that is unoccupied (LUMO)
        # Once we find the first unoccupied orbital, that's the LUMO
        # because the list is sorted by energy.
        if mo.electron_count == 0:
            if lumo is None: # Only assign the first (lowest energy) unoccupied MO
                lumo = mo
            # Optimization: If we found a LUMO and the current MO energy is higher,
            # we don't need to check further for LUMO (though we continue for HOMO).
            # However, simpler to just let the loop finish.

    # Safety check: If LUMO energy is somehow lower than HOMO, something is wrong (e.g., list not sorted)
    if homo and lumo and lumo.energy < homo.energy:
        print(f"Warning: Potential issue finding HOMO/LUMO. LUMO energy ({lumo.energy:.3f}) < HOMO energy ({homo.energy:.3f}). Check MO sorting.", file=sys.stderr)

    return homo, lumo


def calculate_homo_lumo_gap(homo: Optional[MolecularOrbital], lumo: Optional[MolecularOrbital]) -> Optional[float]:
    """
    Calculates the energy gap between HOMO and LUMO.

    Args:
        homo: The HOMO MolecularOrbital object (or None).
        lumo: The LUMO MolecularOrbital object (or None).

    Returns:
        The energy gap in eV if both HOMO and LUMO exist, otherwise None.
    """
    if homo is not None and lumo is not None:
        # Ensure energies are numbers before subtracting
        if isinstance(homo.energy, (int, float)) and isinstance(lumo.energy, (int, float)):
             gap = lumo.energy - homo.energy
             if gap < 0:
                  # This shouldn't happen if find_homo_lumo worked correctly with sorted MOs
                  print(f"Warning: Calculated HOMO-LUMO gap is negative ({gap:.3f} eV). LUMO energy might be lower than HOMO.", file=sys.stderr)
                  return None # Indicate an issue
             return gap
        else:
             print("Warning: HOMO or LUMO energy is not a valid number.", file=sys.stderr)
             return None
    else:
        # Gap is not applicable if either HOMO or LUMO wasn't found
        return None

# --- END NEW FUNCTIONS ---