# mo_diagram_generator/atoms.py
import periodictable
import sys
from typing import Dict, List, Optional, Tuple

# --- Approximate Valence Shell Ionization Potentials (VSIPs) in eV ---
# (Keep VSIP_ENERGIES_EV dictionary as it is - No changes needed here)
VSIP_ENERGIES_EV: Dict[str, float] = {
    # Period 1
    'H_1s': -13.6,
    'He_1s': -24.6, # Noble gas, unlikely to bond simply
    # Period 2
    'Li_2s': -5.4,
    'Be_2s': -9.3,
    'B_2s': -14.0, 'B_2p': -8.3,
    'C_2s': -19.4, 'C_2p': -10.7,
    'N_2s': -25.6, 'N_2p': -13.2,
    'O_2s': -32.3, 'O_2p': -15.9,
    'F_2s': -40.2, 'F_2p': -18.7,
    'Ne_2s': -48.5, 'Ne_2p': -21.6, # Noble gas
    # Period 3 (Example - extend as needed)
    'Na_3s': -5.1,
    'Mg_3s': -7.6,
    'Al_3s': -11.3, 'Al_3p': -6.0,
    'Si_3s': -15.0, 'Si_3p': -7.8,
    'P_3s': -18.7, 'P_3p': -10.0,
    'S_3s': -20.7, 'S_3p': -11.6,
    'Cl_3s': -25.3, 'Cl_3p': -13.7,
    'Ar_3s': -29.2, 'Ar_3p': -15.8, # Noble gas
    # Period 4 Examples (Extend and Refine Energies!)
    'K_4s': -4.3,
    'Ca_4s': -6.1,
    'Sc_4s': -6.5, 'Sc_3d': -7.9, # Approximate
    'Ti_4s': -6.8, 'Ti_3d': -8.9, # Approximate
    'V_4s': -6.7,  'V_3d': -9.0,  # Approximate
    'Cr_4s': -6.8, 'Cr_3d': -10.5, # Updated Cr_3d
    'Mn_4s': -7.4, 'Mn_3d': -11.0, # Updated Mn_3d
    'Fe_4s': -7.9, 'Fe_3d': -12.0, # Updated Fe_3d
    'Co_4s': -7.9, 'Co_3d': -11.7, # Approximate
    'Ni_4s': -7.6, 'Ni_3d': -12.0, # Approximate
    'Cu_4s': -7.7, 'Cu_3d': -12.2, # Approximate (Very Rough!)
    'Zn_4s': -9.4, 'Zn_3d': -15.8, # Approximate
    'Ga_4s': -12.6, 'Ga_4p': -6.0,
    'Ge_4s': -15.6, 'Ge_4p': -7.6,
    'Br_4s': -27.5, 'Br_4p': -14.5, # Added Br (Approx Br_4s needed too)
    'Kr_4s': -32.0, 'Kr_4p': -16.0, # Added Kr (Approx Kr_4s needed too)
}

# Rough energy ordering inversion point for sigma(2p) vs pi(2p)
# Applies typically to B2, C2, N2
SIGMA_PI_INVERSION_Z = 7 # Nitrogen

# Known data for fallback - PRIORITY SOURCE for valence config/electrons
# (Keep _KNOWN_DATA dictionary as it is)
_KNOWN_DATA = {
    # Z: {'period': P, 'group': G, 'valence_e': VE, 'config': ['ns', 'np', ...]}
    1:  {'period': 1, 'group': 1,  'valence_e': 1, 'config': ['1s']},  # H
    2:  {'period': 1, 'group': 18, 'valence_e': 2, 'config': ['1s']},  # He
    3:  {'period': 2, 'group': 1,  'valence_e': 1, 'config': ['2s']},  # Li
    4:  {'period': 2, 'group': 2,  'valence_e': 2, 'config': ['2s']},  # Be
    5:  {'period': 2, 'group': 13, 'valence_e': 3, 'config': ['2s', '2p']}, # B
    6:  {'period': 2, 'group': 14, 'valence_e': 4, 'config': ['2s', '2p']}, # C
    7:  {'period': 2, 'group': 15, 'valence_e': 5, 'config': ['2s', '2p']}, # N
    8:  {'period': 2, 'group': 16, 'valence_e': 6, 'config': ['2s', '2p']}, # O
    9:  {'period': 2, 'group': 17, 'valence_e': 7, 'config': ['2s', '2p']}, # F
    10: {'period': 2, 'group': 18, 'valence_e': 8, 'config': ['2s', '2p']}, # Ne
    11: {'period': 3, 'group': 1,  'valence_e': 1, 'config': ['3s']},  # Na
    12: {'period': 3, 'group': 2,  'valence_e': 2, 'config': ['3s']},  # Mg
    13: {'period': 3, 'group': 13, 'valence_e': 3, 'config': ['3s', '3p']}, # Al
    14: {'period': 3, 'group': 14, 'valence_e': 4, 'config': ['3s', '3p']}, # Si
    15: {'period': 3, 'group': 15, 'valence_e': 5, 'config': ['3s', '3p']}, # P
    16: {'period': 3, 'group': 16, 'valence_e': 6, 'config': ['3s', '3p']}, # S
    17: {'period': 3, 'group': 17, 'valence_e': 7, 'config': ['3s', '3p']}, # Cl
    18: {'period': 3, 'group': 18, 'valence_e': 8, 'config': ['3s', '3p']}, # Ar
    19: {'period': 4, 'group': 1,  'valence_e': 1, 'config': ['4s']},  # K
    20: {'period': 4, 'group': 2,  'valence_e': 2, 'config': ['4s']},  # Ca
    # Transition Metals
    21: {'period': 4, 'group': 3,  'valence_e': 3, 'config': ['4s', '3d']}, # Sc
    22: {'period': 4, 'group': 4,  'valence_e': 4, 'config': ['4s', '3d']}, # Ti
    23: {'period': 4, 'group': 5,  'valence_e': 5, 'config': ['4s', '3d']}, # V
    24: {'period': 4, 'group': 6,  'valence_e': 6, 'config': ['4s', '3d']}, # Cr
    25: {'period': 4, 'group': 7,  'valence_e': 7, 'config': ['4s', '3d']}, # Mn
    26: {'period': 4, 'group': 8,  'valence_e': 8, 'config': ['4s', '3d']}, # Fe
    27: {'period': 4, 'group': 9,  'valence_e': 9, 'config': ['4s', '3d']}, # Co
    28: {'period': 4, 'group': 10, 'valence_e': 10, 'config': ['4s', '3d']},# Ni
    29: {'period': 4, 'group': 11, 'valence_e': 11, 'config': ['4s', '3d']},# Cu
    30: {'period': 4, 'group': 12, 'valence_e': 12, 'config': ['4s', '3d']},# Zn
    # P-block Period 4
    31: {'period': 4, 'group': 13, 'valence_e': 3, 'config': ['4s', '4p']}, # Ga
    32: {'period': 4, 'group': 14, 'valence_e': 4, 'config': ['4s', '4p']}, # Ge
    33: {'period': 4, 'group': 15, 'valence_e': 5, 'config': ['4s', '4p']}, # As
    34: {'period': 4, 'group': 16, 'valence_e': 6, 'config': ['4s', '4p']}, # Se
    35: {'period': 4, 'group': 17, 'valence_e': 7, 'config': ['4s', '4p']}, # Br
    36: {'period': 4, 'group': 18, 'valence_e': 8, 'config': ['4s', '4p']}, # Kr
}


def get_ao_energy(element_symbol: str, orbital_symbol: str) -> Optional[float]:
    """
    Retrieves the approximate VSIP energy for a given atomic orbital.
    (No changes needed in this function's logic)
    """
    key = f"{element_symbol}_{orbital_symbol}"
    energy = VSIP_ENERGIES_EV.get(key)
    if energy is None:
        raise ValueError(f"Energy data for {element_symbol} {orbital_symbol} not found. Update VSIP_ENERGIES_EV in atoms.py.")
    return energy

def get_atom_data(symbol: str) -> Dict:
    """
    Retrieves basic data and valence configuration for an atom.
    (No changes needed in this function's logic, relies on _KNOWN_DATA)
    """
    element = None
    atomic_num = None
    period = None
    group = None
    valence_config_types = None

    try:
        element = periodictable.elements.symbol(symbol)
        atomic_num = element.number
    except (AttributeError, ValueError, KeyError):
        raise ValueError(f"Could not retrieve basic data for symbol '{symbol}' using periodictable library.")

    known_data = _KNOWN_DATA.get(atomic_num, {})

    period = getattr(element, 'period', None) or known_data.get('period')
    if period is None:
        if 1 <= atomic_num <= 2: period = 1
        elif 3 <= atomic_num <= 10: period = 2
        elif 11 <= atomic_num <= 18: period = 3
        elif 19 <= atomic_num <= 36: period = 4
        else: raise ValueError(f"Cannot determine period for Z={atomic_num}.")
        print(f"Warning: Determined period {period} for Z={atomic_num} via fallback range.", file=sys.stderr)

    group = getattr(element, 'group', None) or known_data.get('group')
    if group == 0: group = 18
    if group is None:
         print(f"Warning: Could not determine group for Z={atomic_num}.", file=sys.stderr)

    valence_config_types = known_data.get('config')
    if not valence_config_types:
        print(f"Warning: Valence config not in _KNOWN_DATA for Z={atomic_num}. Inferring PRINCIPAL SHELL {period} only.", file=sys.stderr)
        if period == 1: valence_config_types = [f'{period}s']
        elif period <= 3: valence_config_types = [f'{period}s', f'{period}p']
        elif period == 4:
            print(f"Warning: Period {period} fallback assumes {period}s, {period}p, and possibly {(period-1)}d orbitals. May be incorrect.", file=sys.stderr)
            base_config = [f'{period}s', f'{period}p']
            if group and 3 <= group <= 12:
                 base_config.append(f'{period-1}d')
            valence_config_types = base_config
        else: raise ValueError(f"Valence configuration inference fallback not implemented for period {period}.")

    def orbital_sort_key(s: str) -> Tuple[int, int]:
        n = int(s[0])
        l_char = s[1]
        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        l = l_map.get(l_char, 99)
        return (n, l)

    valence_config_types.sort(key=orbital_sort_key)

    for orb_type in valence_config_types:
        try:
             get_ao_energy(symbol, orb_type)
        except ValueError as e:
             raise ValueError(f"Missing required energy data for atom {symbol}: {e}")

    return {
        "symbol": symbol,
        "atomic_number": atomic_num,
        "period": period,
        "group": group,
        "valence_orbital_types": valence_config_types,
    }


# --- NEW FUNCTION to get neutral valence electron count ---
def get_neutral_valence_electrons(atom_data: Dict) -> int:
    """
    Determines the number of valence electrons for a NEUTRAL atom.
    Uses _KNOWN_DATA first, then falls back to group number.

    Args:
        atom_data: Data dictionary for the atom (from get_atom_data).

    Returns:
        The number of valence electrons for the neutral atom.

    Raises:
        ValueError: If valence electrons cannot be determined.
    """
    atomic_num = atom_data['atomic_number']
    symbol = atom_data['symbol']
    known_data = _KNOWN_DATA.get(atomic_num)

    ve_neutral = -1 # Sentinel value

    # Priority 1: Use explicit valence_e from _KNOWN_DATA
    if known_data and 'valence_e' in known_data:
        ve_neutral = known_data['valence_e']
    else:
        # Priority 2: Fallback using group number
        group = atom_data['group']
        period = atom_data['period']
        if group is None:
            raise ValueError(f"Cannot determine valence electrons for {symbol} (Z={atomic_num}): Missing group info and not in _KNOWN_DATA.")

        if 1 <= group <= 2: # s-block
            ve_neutral = group
        elif 13 <= group <= 18: # p-block
            ve_neutral = group - 10
        elif 3 <= group <= 12: # d-block (Transition Metals)
             print(f"Warning: Using group number ({group}) as neutral valence electron count for {symbol} (Z={atomic_num}, Period {period}). "
                   f"This is often an approximation for d-block elements. Define in _KNOWN_DATA for accuracy.", file=sys.stderr)
             ve_neutral = group
        else:
            raise ValueError(f"Cannot determine valence electrons for {symbol} (Z={atomic_num}) with group {group}. Check data.")

        if ve_neutral != -1:
            print(f"Note: Inferred {ve_neutral} valence electrons for neutral {symbol} using group {group} fallback.", file=sys.stderr)

    if ve_neutral == -1:
         raise ValueError(f"Failed to determine neutral valence electron count for {symbol} (Z={atomic_num}).")

    return ve_neutral


def get_total_valence_electrons(atom_A_data: Dict, atom_B_data: Dict, charge: int) -> int:
    """
    Calculates the total number of valence electrons for the diatomic species,
    ACCOUNTING FOR CHARGE. Uses get_neutral_valence_electrons internally.

    Args:
        atom_A_data: Data dictionary for atom A.
        atom_B_data: Data dictionary for atom B.
        charge: The overall charge of the molecule.

    Returns:
        The total number of valence electrons for the species (atom A + atom B - charge).

    Raises:
        ValueError: If valence electrons cannot be determined or result is negative.
    """
    ve_A_neutral = get_neutral_valence_electrons(atom_A_data)
    ve_B_neutral = get_neutral_valence_electrons(atom_B_data)

    total_ve = ve_A_neutral + ve_B_neutral - charge

    if total_ve < 0:
        raise ValueError(f"Calculated total valence electrons ({total_ve}) is negative. Check input atoms ({atom_A_data['symbol']}, {atom_B_data['symbol']}) and charge ({charge:+}).")

    return total_ve