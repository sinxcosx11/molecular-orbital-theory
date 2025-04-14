# mo_diagram_generator/plotter.py
import matplotlib.pyplot as plt
import numpy as np
from .orbitals import AtomicOrbital, MolecularOrbital
from typing import List, Dict, Tuple, Set, Optional # Added Optional
import math
import sys

# --- Plotting Configuration ---
AO_COLOR = 'cornflowerblue'
MO_COLOR = 'black'
BONDING_COLOR = 'darkblue'
ANTIBONDING_COLOR = 'crimson'
NONBONDING_COLOR = 'dimgrey'
LINE_COLOR = 'darkgrey'
ELECTRON_COLOR = 'black'
AO_LINE_WIDTH = 1.5
MO_LINE_WIDTH = 2.0
CONNECTION_STYLE = ':'
ARROW_STYLE = dict(color=ELECTRON_COLOR, ha='center', va='center', fontsize=12, fontweight='bold')
SUPERSCRIPT_MAP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
UP_ARROW = '↑'
DOWN_ARROW = '↓'

# --- NEW: Helper function to fill AOs for plotting ---
def get_ao_electron_assignments(aos: List[AtomicOrbital], num_valence_electrons: int) -> Dict[int, List[int]]:
    """
    Assigns valence electrons to a list of AOs for plotting purposes.
    Follows Aufbau and Hund's rule roughly.

    Args:
        aos: List of AtomicOrbital objects for ONE atom, assumed sorted by energy.
        num_valence_electrons: Number of valence electrons for that NEUTRAL atom.

    Returns:
        A dictionary mapping the id(ao) to a list of spins [1] or [1, -1].
    """
    assignments = {}
    electrons_remaining = num_valence_electrons

    # Group AOs by energy level for Hund's rule
    levels: Dict[float, List[AtomicOrbital]] = {}
    for ao in aos:
        key = round(ao.energy, 4)
        if key not in levels: levels[key] = []
        levels[key].append(ao)
        assignments[id(ao)] = [] # Initialize empty assignment

    # Fill level by level (Aufbau)
    for energy in sorted(levels.keys()):
        if electrons_remaining <= 0: break
        current_level_aos = levels[energy]
        # Sort within level for consistent filling (e.g., p_0, p_1, p_2)
        current_level_aos.sort(key=lambda ao: ao.symbol)

        # Apply Hund's Rule - Stage 1: Fill singly with spin up
        for ao in current_level_aos:
            if electrons_remaining <= 0: break
            if len(assignments[id(ao)]) == 0: # If empty
                assignments[id(ao)].append(1) # Add spin up
                electrons_remaining -= 1

        if electrons_remaining <= 0: break

        # Apply Hund's Rule - Stage 2: Pair up with spin down
        for ao in current_level_aos:
            if electrons_remaining <= 0: break
            if len(assignments[id(ao)]) == 1: # If singly occupied
                 # Assume capacity 2 for all valence AOs (s, p, d) for simplicity here
                 assignments[id(ao)].append(-1) # Add spin down
                 electrons_remaining -= 1

    if electrons_remaining > 0:
        print(f"Warning (Plotter): Could not assign all {num_valence_electrons} AO electrons (remaining: {electrons_remaining}). Check AO list/valence count.", file=sys.stderr)
    elif electrons_remaining < 0:
        print(f"Warning (Plotter): Overfilled AO electrons (remaining: {electrons_remaining}). Logic error.", file=sys.stderr)


    return assignments


def plot_mo_diagram(
    aos_A: List[AtomicOrbital],
    aos_B: List[AtomicOrbital],
    ve_A: int, # <-- NEW: Valence electrons for neutral Atom A
    ve_B: int, # <-- NEW: Valence electrons for neutral Atom B
    mos: List[MolecularOrbital],
    atom_A_label: str,
    atom_B_label: str,
    title: str,
    bond_order: float,
    magnetic_behavior: str,
    output_filename: str = "mo_diagram.svg"
):
    """Generates and saves an MO diagram using Matplotlib."""

    fig, ax = plt.subplots(figsize=(8, 10)) # Adjusted size slightly if needed
    ax.set_xticks([])
    ax.tick_params(axis='y', direction='in', labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['left'].set_color('darkgrey')
    ax.set_ylabel("Energy (eV)", fontsize=10, labelpad=10)

    # --- Plotting Positions ---
    ao_A_x = 0.15
    ao_B_x = 0.85
    mo_x_center = 0.50
    level_width = 0.08
    # ao_level_width = level_width * 1.2 # Adjust if needed, maybe same as MO width now
    # NEW: Spacing for degenerate AOs (similar to MOs)
    degenerate_ao_slot_spacing = level_width * 1.5
    electron_vertical_offset_factor = 0.015
    electron_x_offset = level_width * 0.35

    # --- Determine Energy Range (includes AOs and MOs) ---
    all_ao_energies = [ao.energy for ao in aos_A + aos_B]
    all_mo_energies = [mo.energy for mo in mos]
    if not all_ao_energies and not all_mo_energies:
        min_energy, max_energy = -20, 10
    else:
        all_energies = all_ao_energies + all_mo_energies
        min_e = min(all_energies) if all_energies else -10
        max_e = max(all_energies) if all_energies else 10
        padding = max(2.0, abs(max_e - min_e) * 0.12)
        min_energy = min_e - padding
        max_energy = max_e + padding * 1.8

    ax.set_ylim(min_energy, max_energy)
    energy_range = max_energy - min_energy
    if energy_range == 0: energy_range = 1.0
    electron_vertical_offset_abs = electron_vertical_offset_factor * energy_range

    # --- Helper to Plot Electrons (no change needed) ---
    def plot_electrons(x_center, y_base, electrons):
        if not electrons: return
        # Ensure electrons is a list, not None
        electrons_list = electrons if electrons is not None else []
        y_pos = y_base + electron_vertical_offset_abs
        if len(electrons_list) == 1:
            arrow = UP_ARROW if electrons_list[0] == 1 else DOWN_ARROW # Should primarily be UP for single
            ax.text(x_center, y_pos, arrow, **ARROW_STYLE)
        elif len(electrons_list) >= 2: # Handle pairs
            # Assume standard pair is [1, -1] after potential sorting
            ax.text(x_center - electron_x_offset, y_pos, UP_ARROW, **ARROW_STYLE)
            ax.text(x_center + electron_x_offset, y_pos, DOWN_ARROW, **ARROW_STYLE)

    # --- Assign AO electrons before plotting ---
    # Sort AOs by energy before assigning electrons
    aos_A.sort(key=lambda ao: ao.energy)
    aos_B.sort(key=lambda ao: ao.energy)
    ao_electrons_A = get_ao_electron_assignments(aos_A, ve_A)
    ao_electrons_B = get_ao_electron_assignments(aos_B, ve_B)

    # --- Plot AOs (MODIFIED FOR DEGENERACY & ELECTRONS) ---
    def plot_ao_levels(x_pos, aos_list, ao_electron_map, label_pos='left'):
        plotted_groups_for_labels: Set[float] = set()
        ao_levels: Dict[float, List[AtomicOrbital]] = {}
        for ao in aos_list:
            key = round(ao.energy, 4)
            if key not in ao_levels: ao_levels[key] = []
            ao_levels[key].append(ao)

        for energy_key, aos_at_level in sorted(ao_levels.items()):
            energy = aos_at_level[0].energy
            # Sort within level for consistent plotting (e.g., p_0, p_1, p_2)
            aos_at_level.sort(key=lambda a: a.symbol)
            num_orbitals = len(aos_at_level)
            total_width = (num_orbitals - 1) * degenerate_ao_slot_spacing if num_orbitals > 1 else 0
            start_x = x_pos - total_width / 2

            # Plot individual AO slots and electrons
            for slot_index, ao in enumerate(aos_at_level):
                slot_center_x = start_x + slot_index * degenerate_ao_slot_spacing
                line_start = slot_center_x - level_width / 2
                line_end = slot_center_x + level_width / 2
                ax.hlines(energy, line_start, line_end, color=AO_COLOR, lw=AO_LINE_WIDTH)

                # Plot AO electrons using the map
                electrons_to_plot = ao_electron_map.get(id(ao), [])
                plot_electrons(slot_center_x, energy, electrons_to_plot)

            # Add Level Labels (group label)
            if energy_key not in plotted_groups_for_labels:
                # Use the base symbol (e.g., "2s", "2p") for the label
                base_symbol = aos_at_level[0].symbol.split('_')[0] # Takes "2p" from "2p_0"
                # Calculate label position relative to the group center
                group_center_x = (start_x + (start_x + total_width)) / 2 # Center of all slots
                text_x = group_center_x + level_width / 2 + 0.018 if label_pos == 'right' else group_center_x - level_width / 2 - 0.018
                ha = 'left' if label_pos == 'right' else 'right'
                ax.text(text_x, energy, base_symbol, ha=ha, va='center', color=AO_COLOR, fontsize=8.5)
                plotted_groups_for_labels.add(energy_key)

    print("Plotting Atom A Orbitals...")
    plot_ao_levels(ao_A_x, aos_A, ao_electrons_A, label_pos='left')
    print("Plotting Atom B Orbitals...")
    plot_ao_levels(ao_B_x, aos_B, ao_electrons_B, label_pos='right')

    # --- Plot MOs (MODIFIED - Ensure Electrons are handled correctly) ---
    # (Logic for MO plotting largely remains the same, but relies on correct electron filling)
    mo_levels: Dict[float, List[MolecularOrbital]] = {}
    for mo in mos:
        key = round(mo.energy, 4)
        if key not in mo_levels: mo_levels[key] = []
        mo_levels[key].append(mo)

    plotted_connections = set()
    processed_mo_groups_for_labels: Set[Tuple[float, str]] = set()

    for energy_key, mos_at_level in sorted(mo_levels.items()):
        energy = mos_at_level[0].energy
        mos_at_level.sort(key=lambda m: m.symbol)
        num_orbitals = len(mos_at_level)
        total_width = (num_orbitals - 1) * degenerate_ao_slot_spacing if num_orbitals > 1 else 0 # Use same spacing? maybe different var needed
        mo_degenerate_slot_spacing = level_width * 1.5 # Define MO specific spacing
        total_width = (num_orbitals - 1) * mo_degenerate_slot_spacing if num_orbitals > 1 else 0
        start_x = mo_x_center - total_width / 2

        # Plot MO slots and electrons
        for slot_index, mo in enumerate(mos_at_level):
            slot_center_x = start_x + slot_index * mo_degenerate_slot_spacing
            mo_line_color = MO_COLOR
            if mo.is_bonding is True: mo_line_color = BONDING_COLOR
            elif mo.is_bonding is False: mo_line_color = ANTIBONDING_COLOR
            elif mo.is_nonbonding: mo_line_color = NONBONDING_COLOR

            ax.hlines(energy, slot_center_x - level_width / 2, slot_center_x + level_width / 2, color=mo_line_color, lw=MO_LINE_WIDTH)
            # Plot MO electrons (ensure mo.electrons list is correctly populated by filler)
            plot_electrons(slot_center_x, energy, mo.electrons)

            # Draw Connecting Lines (Connect to original AO positions)
            for origin_ao in mo.origin_AOs:
                # Determine the X position of the originating AO level group
                ao_x_base = ao_A_x if origin_ao.atom_label == 'A' else ao_B_x
                ao_energy = origin_ao.energy
                # Find the group of AOs at this energy for the specific atom
                ao_group = [ao for ao in (aos_A if origin_ao.atom_label == 'A' else aos_B) if math.isclose(ao.energy, ao_energy, rel_tol=1e-5)]
                num_ao_deg = len(ao_group)
                ao_total_width = (num_ao_deg - 1) * degenerate_ao_slot_spacing if num_ao_deg > 1 else 0
                ao_start_x = ao_x_base - ao_total_width / 2
                # Find the specific slot index of the origin_ao within its degenerate group
                ao_group.sort(key=lambda a: a.symbol)
                try:
                     ao_slot_index = ao_group.index(origin_ao)
                     ao_slot_center_x = ao_start_x + ao_slot_index * degenerate_ao_slot_spacing
                except ValueError:
                     # Fallback if exact AO object not found (shouldn't happen ideally)
                     print(f"Warning: Origin AO {origin_ao.symbol} not found in plotted group. Connecting to base AO X.", file=sys.stderr)
                     ao_slot_center_x = ao_x_base # Connect to the central position

                # Use unique keys including object IDs or precise coords/symbols if needed
                ao_key = (round(ao_energy, 3), round(ao_slot_center_x, 3), id(origin_ao))
                mo_key = (round(energy, 3), round(slot_center_x, 3), id(mo))

                conn_key = tuple(sorted((ao_key, mo_key)))
                if conn_key not in plotted_connections:
                    ax.plot([ao_slot_center_x, slot_center_x], [ao_energy, energy], ls=CONNECTION_STYLE, color=LINE_COLOR, lw=0.8, zorder=-1)
                    plotted_connections.add(conn_key)

        # Add MO Level Labels (handle degeneracy - logic unchanged)
        group_representatives: Dict[str, MolecularOrbital] = {}
        group_electrons: Dict[str, int] = {}
        group_x_coords: Dict[str, List[float]] = {}

        for slot_index, mo in enumerate(mos_at_level):
            group_id = mo.degeneracy_group if mo.degeneracy_group else mo.symbol
            slot_center_x = start_x + slot_index * mo_degenerate_slot_spacing
            if group_id not in group_representatives:
                group_representatives[group_id] = mo
                group_electrons[group_id] = 0
                group_x_coords[group_id] = []
            group_electrons[group_id] += mo.electron_count
            group_x_coords[group_id].append(slot_center_x)

        label_x_positions = []
        labels_to_render = []
        for group_id, rep_mo in group_representatives.items():
            label_key = (energy_key, group_id)
            if label_key not in processed_mo_groups_for_labels:
                coords = group_x_coords[group_id]
                label_x = (min(coords) + max(coords)) / 2 + level_width / 2 + 0.018
                label_x_positions.append(label_x)
                total_electrons = group_electrons[group_id]
                base_symbol = rep_mo.symbol
                if rep_mo.degeneracy_group: base_symbol = base_symbol.rsplit('_', 1)[0]
                electron_str = str(total_electrons).translate(SUPERSCRIPT_MAP) if total_electrons > 0 else ""
                full_label = f"{base_symbol}{electron_str}"
                labels_to_render.append({'x': label_x, 'text': full_label})
                processed_mo_groups_for_labels.add(label_key)

        max_label_x = max(label_x_positions) if label_x_positions else mo_x_center + level_width / 2 + 0.018
        for lbl in labels_to_render:
             ax.text(max_label_x, energy, lbl['text'], ha='left', va='center', color=MO_COLOR, fontsize=8.5)


    # --- Add Atom Labels and Title (remains the same) ---
    label_y_pos = max_energy - padding * 0.15
    ax.text(ao_A_x, label_y_pos, f"Atom A\n{atom_A_label}", ha='center', va='top', fontsize=10, fontweight='medium')
    ax.text(ao_B_x, label_y_pos, f"Atom B\n{atom_B_label}", ha='center', va='top', fontsize=10, fontweight='medium')
    ax.text(mo_x_center, label_y_pos, "MOs", ha='center', va='top', fontsize=10, fontweight='medium')

    bond_order_str = f"{bond_order:.1f}".rstrip('0').rstrip('.') if '.' in f"{bond_order:.1f}" else str(int(bond_order))
    formatted_behavior = magnetic_behavior.replace('e-', 'e⁻')
    full_title = f"MO Diagram for {title}\nBond Order = {bond_order_str} | {formatted_behavior}"
    ax.set_title(full_title, fontsize=14, pad=30, fontweight='bold')

    # --- Final Touches & Save (remains the same) ---
    plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    try:
        output_format = 'svg'
        if '.' in output_filename:
             ext = output_filename.split('.')[-1].lower()
             if ext in ['png', 'svg', 'pdf', 'jpg', 'jpeg']: output_format = ext
             else: output_filename = f"{output_filename.rsplit('.', 1)[0]}.svg"
        else: output_filename = f"{output_filename}.svg"
        plt.savefig(output_filename, bbox_inches='tight', dpi=150, format=output_format)
        print(f"Diagram saved to '{output_filename}'")
    except Exception as e:
        print(f"Error saving plot to '{output_filename}': {e}", file=sys.stderr)
    finally:
        plt.close(fig)