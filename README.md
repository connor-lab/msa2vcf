# msa2vcf

Turn a fasta-format MSA into one vcf per sequence.

`msa2vcf msa.fasta ref_name`

VCF position is based on reference position.

IUPAC bases:
  - Ref base is stripped from list of potentials
  - Variant AF is given as 1/IUPAC redundancy
  
 Pretty much all Python stdlib and all that entails - no input checking, no brakes.
