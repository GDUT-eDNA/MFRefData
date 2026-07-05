import os
import argparse
from pathlib import Path
from db import download_NCBI
from remove_primers import isPCR_optimized
from cluster import cluster
from totax import totax
from filter import filter
from output import output
from vistual_bar import bar
from vistual_pie import pie
from vistual_venn import venn


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Download, process, and generate reference sequence datasets."
    )

    parser.add_argument(
        "-list",
        "--list",
        dest="species_list",
        required=True,
        help=(
            "Species list file. A bare filename such as 'example2' "
            "will be interpreted as './work/example2.txt'."
        ),
    )

    parser.add_argument(
        "-region",
        "--region",
        dest="species_region",
        required=True,
        choices=["12s", "18s", "16s", "COI"],
        help="Target gene region."
    )

    parser.add_argument(
        "-type",
        "--type",
        dest="output_type",
        required=True,
        choices=["qiime2"],
        help="Output data format."
    )

    return parser.parse_args()


args = parse_arguments()
work_dir = './work'
Path(work_dir).mkdir(parents=True, exist_ok=True)

Species_list = args.species_list
Species_region = args.species_region
Output_type = args.output_type

download_NCBI(
    path=os.path.join(work_dir, Species_list),
    region=Species_region,
    count=100
)
isPCR_optimized(
    input_path=os.path.join(work_dir, 'Download'),
    F='GTCGGTAAAACTCGTGCCAGC',
    R='CATAGTGGGGTATCTAATCCCAGTTTG',
    mis=1
)

cluster(
    file_path=os.path.join(work_dir, 'ISPCR', 'ISPCR_sequence'),
    csv_path=os.path.join(work_dir, 'ISPCR', 'ISPCR_info.csv'),
    threshold=0.97
)
totax(
    csv=os.path.join(work_dir, 'Cluster', 'Cluster_info.csv'),
    dump=r'D:\py_problem\Apocalypse\MFRefData\Project\taxdmp'
)
filter(
    csv_file=os.path.join(work_dir, 'Totax', 'Totax.csv'),
    fasta_file=os.path.join(work_dir, 'Cluster', 'Combined_Cluster.fasta'),
    raw_max_len=5000
)
output(
    type=Output_type,  # qiime2 sintax RDP BLAST kraken2 DADA2 mothur
    file=os.path.join(work_dir, 'Filter_5000', 'Filter.fasta'),
    refcsv=os.path.join(work_dir, 'Filter_5000', 'Filter.csv')
)
os.makedirs(os.path.join(work_dir, 'output', 'output_bar'), exist_ok=True)
os.makedirs(os.path.join(work_dir, 'output', 'output_pie'), exist_ok=True)
bar(inutfile=os.path.join(work_dir, 'Filter_5000', 'Filter.csv'),
    output_file=os.path.join(work_dir, 'output', 'output_bar', 'bar.png'))

pie(file=os.path.join(work_dir, 'Filter_5000', 'Filter.csv'),
    output_fp_csv=os.path.join(work_dir, 'output', 'output_pie', 'species_condition.csv'),
    output_fp_fig=os.path.join(work_dir, 'output', 'output_pie', 'pie.png'))
