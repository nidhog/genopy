record_file = "input/dna.example.fasta"
default_record_file_type = "fasta"
start_codons = ['ATG']
end_codons = ['TAA', 'TAG', 'TGA']


def is_start_codon(codon):
    return codon in start_codons

def is_end_codon(codon):
    return codon in end_codons
