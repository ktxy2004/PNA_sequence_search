import multiprocessing as mp
from functools import partial
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import reverse_complement


def find_sequence(genome_file, sequence, gff_file,num_processes=10):
    # Load genome sequences
    genome_records = list(SeqIO.parse(genome_file, "fasta"))

    # Create a dictionary of CDS and exon locations keyed by start and end positions
    cds_dict = {}
    with open(gff_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                feature = cols[2]
                if feature == "CDS" or feature == "exon":
                    start, end = map(int, cols[3:5])
                    strand = cols[6]
                    gene_name = cols[8].split(";")[2]+ "_" + cols[8].split(";")[3] #gene_name = cols[8].split(";")[0].split("=")[1]
                    if feature != "region":
                        cds_dict[(start, end, cols[0])] = (strand, gene_name)

    # Find the locations of the query sequence
    query_length = len(sequence)
    locations = []
    for genome_record in genome_records:
        genome_seq = genome_record.seq
        chromosome_name = genome_record.id
        for strand, nucleotide_seq in [(+1, sequence), (-1, reverse_complement(sequence))]:
            for i in range(len(genome_seq) - query_length + 1):
                num_matching_bases = sum([genome_seq[i+j].upper() == nucleotide_seq[j].upper() for j in range(query_length)])
                if num_matching_bases == query_length:
                    start_pos = i+1
                    end_pos = i+query_length
                    locations.append((chromosome_name, start_pos, end_pos, strand, nucleotide_seq))

    # Find the CDSs that overlap with the query sequence
    overlapping_cds = []
    for cds_start, cds_end, cds_chromosome in cds_dict.keys():
        for query_chromosome, query_start, query_end, query_strand, query_seq in locations:
            if cds_start <= query_end and cds_end >= query_start and query_chromosome == cds_chromosome and cds_dict[(cds_start, cds_end, cds_chromosome)][1] != "region":
                if query_strand == 1:
                    if cds_dict[(cds_start, cds_end, cds_chromosome)][0] == "+":
                        strand_comment = "plus plus match"
                    else:
                        strand_comment = "reverse complement match"
                else:
                    if cds_dict[(cds_start, cds_end, cds_chromosome)][0] == "+":
                        strand_comment = "reverse complement match"
                    else:
                        strand_comment = "plus plus match"
                cds_comment = ""
                if query_start <= cds_start <= query_end or query_start <= cds_end <= query_end:
                    cds_comment = "Translational Start Site Included"
                overlapping_cds.append((cds_dict[(cds_start, cds_end, cds_chromosome)][1], cds_start, cds_end, cds_dict[(cds_start, cds_end, cds_chromosome)][0], query_chromosome, query_start, query_end, query_strand, strand_comment, cds_comment, query_seq))


            
    return overlapping_cds


# Set the filenames and query sequence
#genome_file = "/home/kathrine/Apisum_Buchnera_genome_gff/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.fna" ##AL4f
#gff_file = "/home/kathrine/Apisum_Buchnera_genome_gff/GCF_005508785.2_pea_aphid_22Mar2018_4r6ur_v2_genomic.gff"
#sequence = "GCCATTTGAC"

# Set the filenames and query sequence
genome_file = "/home/kathrine/Apisum_Buchnera_genome_gff/GCF_000142985.2_Acyr_2.0_genomic.fna" ##LSR1
gff_file = "/home/kathrine/Apisum_Buchnera_genome_gff/genomic.gff"
sequence = "GCCATTTGAC"

# Set the filenames and query sequence
#genome_file = "/home/kathrine/Apisum_Buchnera_genome_gff/GCF_000009605.1_ASM960v1_genomic.fna" ##Buchnera
#gff_file = "/home/kathrine/Apisum_Buchnera_genome_gff/GCF_000009605.1_ASM960v1_genomic.gff"
#sequence = "GCCATTTGAC"


def write_results_to_file(results, output_file):
    with open(output_file, "w") as f:
        f.write("Gene Name\tCDS start position\tCDS end position\tStrand\tQuery start position\tQuery end position\tQuery strand\tComment\tExtra Comment\tQuery sequence\n")
        for result in results:
            f.write("\t".join(map(str, result)) + "\n")

# Run the function and write the results to a tab-delimited file
output_file = "/home/kathrine/Apisum_Buchnera_genome_gff/20230411_bias5_LSR1_GCCATTTGAC_output.txt"
results = find_sequence(genome_file, sequence, gff_file)
write_results_to_file(results, output_file)

