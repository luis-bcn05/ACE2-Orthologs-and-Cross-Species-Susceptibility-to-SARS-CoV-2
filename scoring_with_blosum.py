from Bio.Align import substitution_matrices

def read_fasta(file_path):
    sequences = {}
    current_id = None
    current_seq = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                
                # Guardar la secuencia anterior si existe
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                # Extraer ID o nombre
                current_id = line[1:].split()[0]  # e.g., "NP_001358344.1"
                current_seq = []
            else:
                current_seq.append(line)
        # Guardar la última secuencia
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences

sequences = read_fasta("ACE2_refseq_protein.fasta")

# Load full BLOSUM62 matrix
blosum62 = substitution_matrices.load("BLOSUM62")

# Your dictionary of positions and expected amino acids (1-based positions)
expected_residues = {
    23: 'Q',
    29: 'D',
    30: 'K',
    33: 'H',
    34: 'E',
    40: 'Y',
    41: 'Q',
    44: 'L',
    81: 'M',
    82: 'Y',
    352: 'K',
    353: 'G',
    354: 'D',
    37: 'D',
    392: 'R',
    89: 'N',
}

def score_sequences_blosum_positions(sequences, expected_residues):
    scores = {}
    for seq_id, seq in sequences.items():
        score = 0
        for pos_1based, expected_aa in expected_residues.items():
            pos0 = pos_1based   # convert to zero-based indexing for Python
            if pos0 >= len(seq):
                # Position missing, you can assign penalty or skip
                # For example, penalize by -4 (like a strong mismatch)
                score += -4
                continue
            observed_aa = seq[pos0]
            if observed_aa == '-':
                # Gap, also penalize
                score += -4
                continue
            try:
                score += blosum62[expected_aa, observed_aa]
            except KeyError:
                # Unknown amino acid encountered, assign a penalty
                score += -4
        scores[seq_id] = score
    return scores

scores = score_sequences_blosum_positions(sequences, expected_residues)

for seq_id, score in sorted(scores.items(), key=lambda x: x[1], reverse=True):
    print(f"{seq_id}\t{score}")
