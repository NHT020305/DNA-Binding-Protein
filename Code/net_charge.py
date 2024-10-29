import matplotlib.pyplot as plt



####################  NET CHARGE  ####################


selected_pH = 7.0


pKa_values = {
        'K': 10.5,  # Lysine
        'R': 12.5,  # Arginine
        'H': 6.0,   # Histidine
        'D': 3.9,   # Aspartic acid
        'E': 4.2    # Glutamic acid
    }


def net_charge(protein_sequence, pH=selected_pH):
    
    charge = 0
    
    for amino_acid in protein_sequence:
        if amino_acid in pKa_values:
            if amino_acid in ['K', 'R', 'H']:  # Basic amino acids
                if pH < pKa_values[amino_acid]:
                    charge += 1  # Protonated form (positive charge)
            elif amino_acid in ['D', 'E']:  # Acidic amino acids
                if pH > pKa_values[amino_acid]:
                    charge -= 1  # Deprotonated form (negative charge)
    
    return charge

