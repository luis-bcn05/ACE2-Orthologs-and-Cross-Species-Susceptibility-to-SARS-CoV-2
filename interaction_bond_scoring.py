
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
# Define interaction groups and weights
interaction_profiles = {
    'hydrophobic': {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y'},
    'hbond': {'S', 'T', 'N', 'Q', 'Y', 'H', 'E', 'D', 'K', 'R'},
    'polar': {'S', 'T', 'N', 'Q', 'Y', 'C', 'H'},
    'electrostatic': {'K', 'R', 'H', 'D', 'E'},
    'salt': {'D', 'E', 'K', 'R', 'H'},
    'backbone': {'G', 'A'}
}

interaction_weights = {
    'hydrophobic': 3.0,
    'hbond': 2.5,
    'polar': 2.0,
    'electrostatic': 2.0,
    'salt': 1.5,
    'backbone': 1.0,
    'glycan': 1.0
}

# Classify what interaction types a residue supports
def classify_residue(residue):
    supported = set()
    for itype, residues in interaction_profiles.items():
        if residue in residues:
            supported.add(itype)
    print(f"{residue} supports: {supported}")

    return supported

# Score sequences
def score_sequences_weighted(sequences, interaction_map):
    """
    sequences: list of protein sequences (same length, aligned)
    interaction_map: dict {position: [interaction_type1, interaction_type2, ...]}
    
    Returns: list of scores
    """
    scores = []
    for seq in sequences:
        score = 0.0
        for pos, interactions in interaction_map.items():
            if pos >= len(sequences[seq]):
                score += sum(interaction_weights[i] for i in interactions)  # penalize missing
                # print(f"{seq} missing pos {pos}")

                continue
            residue = sequences[seq][pos]
            capabilities = classify_residue(residue)
            for interaction in interactions:
                if interaction not in capabilities:
                    score += interaction_weights.get(interaction, 1.0)
        scores.append(score)
    return scores

# Example usage with real interaction map from the PDF
interaction_map = {
    23: ['hbond', 'hydrophobic'],    # Q24
    29: ['salt', 'polar'],           # D30
    30: ['electrostatic', 'hbond'],  # K31
    33: ['polar', 'hbond'],          # H34
    34: ['hbond'],                   # E35
    40: ['hbond'],                   # Y41
    41: ['polar'],                   # Q42
    44: ['hydrophobic'],             # L45
    81: ['hydrophobic'],             # M82
    82: ['hbond'],                   # Y83
    352: ['hbond', 'electrostatic'], # K353
    353: ['backbone'],              # G354
    354: ['electrostatic'],          # D355
    37: ['hbond'],                   # D38
    392: ['electrostatic', 'polar'], # R393
    89: ['glycan']                   # N90 
}

# Test sequences (use aligned ones)
# sequences = [
#     "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFP" + "A" * 320,  # Add padding to simulate full length
#     "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGFFA" + "A" * 320,
#     "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYLA" + "A" * 320
# ]

scores = score_sequences_weighted(sequences, interaction_map)
# print(scores)
print(sequences)

seq_scores = list(zip(sequences.keys(), scores))

# Sort by score descending
seq_scores_sorted = sorted(seq_scores, key=lambda x: x[1], reverse=True)

with open('results.txt', 'w') as f:
    for seq_id, score in seq_scores_sorted:
        f.write(f"{seq_id}\t{score}\n")


with open('prueba.txt', 'w') as f:
    for seq_id, seq in sequences.items():
        f.write(f"{seq_id}: {seq[:50]}...")  # Print first 50 characters for brevity
# Print sequences for verification

with open('debugging.txt', 'w') as f:
    for sid, seq in sequences.items():
        f.write(f"\n")
        f.write(f"\n{sid}")
        for pos in sorted(interaction_map.keys()):
            if pos < len(seq):
                f.write(f"  Pos {pos+1}: {seq[pos]}")
            else:
                f.write(f"  Pos {pos+1}: [out of range]")



### i point is that i dont have a reference sequence, i have some positions that tells me what is the aminoacid that should go there.