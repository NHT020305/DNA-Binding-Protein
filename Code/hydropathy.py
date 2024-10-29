import matplotlib.pyplot as plt


####################  HYDROPATHY  ####################


hydropathy_index = {
    'A': 1.8,   # Alanine
    'C': 2.5,   # Cysteine
    'D': -3.5,  # Aspartic Acid
    'E': -3.5,  # Glutamic Acid
    'F': 2.8,   # Phenylalanine
    'G': -0.4,  # Glycine
    'H': -0.2,  # Histidine
    'I': 4.5,   # Isoleucine
    'K': -3.9,  # Lysine
    'L': 3.8,   # Leucine
    'M': 1.9,   # Methionine
    'N': -3.5,  # Asparagine
    'P': -1.6,  # Proline
    'Q': -3.5,  # Glutamine
    'R': -4.5,  # Arginine
    'S': -0.8,  # Serine
    'T': -0.7,  # Threonine
    'V': 4.2,   # Valine
    'W': -0.9,  # Tryptophan
    'Y': -1.3,  # Tyrosine
}


def hydropathy(protein_sequence, amino_acid_list):
    hydropathy = 0
    for amino_acid in protein_sequence:
        if amino_acid not in amino_acid_list:
            continue
        hydropathy += hydropathy_index[amino_acid]
    return hydropathy / len(protein_sequence) 
