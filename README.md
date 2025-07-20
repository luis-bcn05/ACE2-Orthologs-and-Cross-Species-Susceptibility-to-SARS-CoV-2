# 🧬 ACE2 Orthologs and Cross-Species Susceptibility to SARS-CoV-2

In this bioinformatics project, we explored the evolutionary diversity of the ACE2 protein across mammalian species to investigate how variations in the sequence could affect susceptibility to SARS-CoV-2 infection.

SARS-CoV-2 uses the ACE2 receptor to enter human cells. Since ACE2 is a membrane-bound protein expressed in many vertebrates, comparing its sequence across species can help us understand why some animals (like cats or ferrets) are susceptible to COVID-19 while others (like dogs or rodents) show resistance.

---

## 🔍 Objectives

- Analyze ACE2 orthologs across placentals, marsupials, and monotremes
- Perform multiple sequence alignment and phylogenetic tree construction
- Extract and compare 16 key residues critical for SARS-CoV-2 binding
- Investigate susceptibility by comparing scoring strategies based on:
  - **Interaction bond contribution**
  - **BLOSUM62 substitution matrix**

---

## 🧪 Methodology

1. **Sequence Retrieval**  
   Protein sequences were collected from the NCBI database (in FASTA format), including an outgroup (chicken) for rooting the phylogenetic tree.

2. **Multiple Sequence Alignment (MSA)**  
   MSA was performed using **MAFFT** with L-INS-i mode due to its high accuracy and ability to handle divergent sequences.

3. **Phylogenetic Tree Construction**  
   Two trees were generated:
   - Based on the full ACE2 protein (805 amino acids)
   - Based on a 16-residue subset involved in SARS-CoV-2 binding (extracted using a custom Python function)

4. **Susceptibility Scoring**  
   Two Python scripts are included in this repository:
   - `score_interaction.py`: Scores sequences based on interaction bond contribution of critical ACE2 residues
   - `score_blosum.py`: Uses the BLOSUM62 matrix to assess similarity of critical residues to human ACE2
  

---

## 🧠 Key Insights

- Variability in ACE2 binding residues correlates with species susceptibility to SARS-CoV-2.
- Relying solely on overall sequence similarity may be misleading; functional residues matter most.
- Other biological factors (e.g., ACE2 expression, TMPRSS2 presence) also contribute to viral entry resistance.

---

## 👥 Authors

- Luis Carlos Ospina  
- David Lázaro  
- Jan Jardí

---

## 📌 Notes

This project was developed during the second year of the Bioinformatics degree as part of a team research and presentation assignment. The analysis combines comparative genomics, evolutionary biology, and Python scripting to offer a functional view of viral-host interactions across species.

---

## 📜 License

This repository is shared for educational and academic purposes. Please credit the authors if reused.


