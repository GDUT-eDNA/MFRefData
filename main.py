import os
from db import download_NCBI
from remove_primers import isPCR_optimized
from cluster import cluster
from totax import totax
from filter import filter
from output import output
from vistual_bar import bar
from vistual_pie import pie
from vistual_venn import venn

work_dir = './work'

yp_str = "7"
if yp_str == "1":
    download_NCBI(
        path=os.path.join(work_dir, '1.txt'),
        region='12s',
        count=100
    )
elif yp_str == "2":
    isPCR_optimized(
        input_path=os.path.join(work_dir, 'Download'),
        F='GTCGGTAAAACTCGTGCCAGC',
        R='CATAGTGGGGTATCTAATCCCAGTTTG',
        mis=1
    )
elif yp_str == "3":
    cluster(
        file_path=os.path.join(work_dir, 'ISPCR', 'ISPCR_sequence'),
        csv_path=os.path.join(work_dir, 'ISPCR', 'ISPCR_info.csv'),
        threshold=0.97
    )
elif yp_str == "4":
    totax(
        csv=os.path.join(work_dir, 'Cluster', 'Cluster_info.csv'),
        dump=r'D:\py_problem\Apocalypse\MFRefData\Project\taxdmp'
    )
    #https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz taxdmp下载链接
elif yp_str == "5":
    filter(
        csv_file=os.path.join(work_dir, 'Totax','Totax.csv'),
        fasta_file=os.path.join(work_dir, 'Cluster', 'Combined_Cluster.fasta'),
        raw_max_len=5000
    )
elif yp_str == "6":
    output(
        type='qiime2', #qiime2 sintax RDP BLAST kraken2 DADA2 mothur
        file=os.path.join(work_dir, 'Filter_5000', 'Filter.fasta'),
        refcsv=os.path.join(work_dir, 'Filter_5000', 'Filter.csv')
    )
elif yp_str == "7":
    os.makedirs(os.path.join(work_dir, 'output','output_bar'), exist_ok=True)
    os.makedirs(os.path.join(work_dir, 'output', 'output_pie'), exist_ok=True)
    os.makedirs(os.path.join(work_dir, 'output', 'output_venn'), exist_ok=True)

    bar(inutfile=os.path.join(work_dir, 'Filter_5000', 'Filter.csv'),
        output_file=os.path.join(work_dir, 'output','output_bar', 'bar.png'))

    pie(file=os.path.join(work_dir, 'Filter_5000', 'Filter.csv'),
        output_fp_csv=os.path.join(work_dir, 'output', 'output_pie', 'species_condition.csv'),
        output_fp_fig=os.path.join(work_dir, 'output', 'output_pie', 'pie.png'))

    venn(
        file=os.path.join(work_dir, 'Filter_5000', 'Filter.csv'),
        species_list=os.path.join(work_dir, '1.txt'),
        formation='qiime2',
        format_path=r'D:\py_problem\Apocalypse\bioinformatics\data_share',
        output_dir=os.path.join(work_dir, 'output', 'output_venn')
    )
else: print("wrong!")


