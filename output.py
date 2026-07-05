import os.path
import os
import re
import copy

from Bio import SeqIO
import pandas as pd
import numpy as np
from pathlib import Path


def remove_prefix(value):
    value = str(value)
    if len(value) > 2 and value[1] == '_':
        return value[2:]
    return value


def read_refcsv(refcsv):
    return pd.read_csv(refcsv, sep=None, engine='python', encoding='utf_8_sig')


def output_qiime(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output qiime2')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa = os.path.join(output_file, 'output_qiime2.fasta')
            output_file_txt = os.path.join(output_file, 'output_qiime2.txt')

            txt_format = []
            fasta_format = []
            csv_df = pd.read_csv(refcsv, encoding='utf_8_sig')
            seq_name = csv_df['name'].value_counts().index.tolist()
            seq_name_with_tax = csv_df[['name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]
            # 先处理序列的新id
            seq_record_dict = {}
            for seq_record in SeqIO.parse(file, 'fasta'):
                # print(seq_record)
                seq_record_dict[seq_record.id] = seq_record
            #
            for name in seq_name:
                if name in seq_record_dict.keys():
                    seq = seq_record_dict[name]
                    fasta_format.append(seq)
                else:
                    print(f'{name} is not in sequence')
                    continue

                # 处理后续信息
                new_tax = seq_name_with_tax[seq_name_with_tax['name'] == name].iloc[0, 1:].tolist()
                new_tax = ';'.join(new_tax)
                new_tax = re.sub('_', repl='__', string=new_tax)

                new_tax_df = f'{name}	{new_tax}'
                txt_format.append(new_tax_df)

            SeqIO.write(fasta_format, output_file_fa, 'fasta')

            with open(output_file_txt, 'w', encoding='utf_8_sig') as f:
                f.write('\n'.join(txt_format))


def output_sintax(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output sintax')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa = os.path.join(output_file, 'output_sintax.fasta')

            df = read_refcsv(refcsv)
            fasta_sintax = []
            # 使用name建立联系
            for seq_record in SeqIO.parse(file, 'fasta'):
                name = seq_record.id
                print(name)
                for index, row in df.iterrows():
                    if row['name'] == name:
                        # 处理序列信息
                        kingdom = remove_prefix(row['kingdom'])
                        phylum = remove_prefix(row['phylum'])
                        class_name = remove_prefix(row['class'])
                        order = remove_prefix(row['order'])
                        family = remove_prefix(row['family'])
                        genus = remove_prefix(row['genus'])
                        species = remove_prefix(row['species']).replace(' ', '_')

                        seq_record.id = f"{row['name']};tax=k:{kingdom},p:{phylum},c:{class_name},o:{order},f:{family},g:{genus},s:{species};"
                        seq_record.description = ''
                        fasta_sintax.append(seq_record)

            SeqIO.write(fasta_sintax, output_file_fa, 'fasta')
            print('finish sintax output')


def output_RDP(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output RDP')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa = os.path.join(output_file, 'output_RDP.fasta')

            df = read_refcsv(refcsv)
            fasta_sintax = []
            # 使用name建立联系
            for seq_record in SeqIO.parse(file, 'fasta'):
                name = seq_record.id
                print(name)
                for index, row in df.iterrows():
                    if row['name'] == name:
                        # 处理序列信息
                        kingdom = remove_prefix(row['kingdom'])
                        phylum = remove_prefix(row['phylum'])
                        class_name = remove_prefix(row['class'])
                        order = remove_prefix(row['order'])
                        family = remove_prefix(row['family'])
                        genus = remove_prefix(row['genus'])

                        seq_record.id = f"{row['name']} Root;{kingdom};{phylum};{class_name};{order};{family};{genus}"
                        seq_record.description = ''
                        fasta_sintax.append(seq_record)

            SeqIO.write(fasta_sintax, output_file_fa, 'fasta')
            print('finish RDP output')


def output_BLAST(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output BLAST')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa = os.path.join(output_file, 'output_BLAST.fasta')
            output_file_tsv = os.path.join(output_file, 'output_BLAST.tsv')

            df = read_refcsv(refcsv)
            new_df = []
            fasta_sintax = []

            # 使用name建立联系
            for seq_record in SeqIO.parse(file, 'fasta'):
                name = seq_record.id
                print(name)
                for index, row in df.iterrows():
                    if row['name'] == name:
                        seq_record.description = ''
                        fasta_sintax.append(seq_record)
                        new_df.append([row['name'], row['taxid']])

            SeqIO.write(fasta_sintax, output_file_fa, 'fasta')
            new_df = pd.DataFrame(new_df)
            new_df.to_csv(output_file_tsv, sep='\t', index=False, header=False, encoding='utf-8')
            print('finish BLAST output')


def output_kraken2(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output kraken2')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa = os.path.join(output_file, 'output_kraken2.fasta')

            df = read_refcsv(refcsv)
            fasta_sintax = []
            # 使用name建立联系
            for seq_record in SeqIO.parse(file, 'fasta'):
                name = seq_record.id
                print(name)
                for index, row in df.iterrows():
                    if row['name'] == name:
                        # 处理序列信息
                        species = remove_prefix(row['species']).replace(' ', '_')

                        seq_record.id = f"{row['name']}|kraken:taxid|{row['taxid']} {species}"
                        seq_record.description = ''
                        fasta_sintax.append(seq_record)

            SeqIO.write(fasta_sintax, output_file_fa, 'fasta')
            print('finish kraken2 output')


def output_DADA2(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output DADA2')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa1 = os.path.join(output_file, 'output_DADA2_full.fasta')
            output_file_fa2 = os.path.join(output_file, 'output_DADA2_genus_species.fasta')

            df = read_refcsv(refcsv)
            fasta_sintax1 = []
            fasta_sintax2 = []
            # 使用name建立联系
            for seq_record in SeqIO.parse(file, 'fasta'):
                name = seq_record.id
                print(name)
                for index, row in df.iterrows():
                    if row['name'] == name:
                        # 处理序列信息
                        kingdom = remove_prefix(row['kingdom'])
                        phylum = remove_prefix(row['phylum'])
                        class_name = remove_prefix(row['class'])
                        order = remove_prefix(row['order'])
                        family = remove_prefix(row['family'])
                        genus = remove_prefix(row['genus'])
                        species = remove_prefix(row['species']).replace('_', ' ')

                        seq_record1 = copy.deepcopy(seq_record)
                        seq_record1.id = f"{kingdom};{phylum};{class_name};{order};{family};{genus};{species};"
                        seq_record1.description = ''
                        fasta_sintax1.append(seq_record1)

                        seq_record2 = copy.deepcopy(seq_record)
                        if species.startswith(genus + ' '):
                            seq_record2.id = f"{row['name']} {species}"
                        else:
                            seq_record2.id = f"{row['name']} {genus} {species}"
                        seq_record2.description = ''
                        fasta_sintax2.append(seq_record2)

            SeqIO.write(fasta_sintax1, output_file_fa1, 'fasta')
            SeqIO.write(fasta_sintax2, output_file_fa2, 'fasta')
            print('finish DADA2 output')


def output_mothur(fasta_file, refcsv):
    file = fasta_file
    refcsv = refcsv
    if os.path.exists(file) == True:
        if os.path.exists(refcsv) == True:
            print('output mothur')
            output_dir = os.path.dirname(os.path.dirname(file))
            output_file = os.path.join(output_dir, 'output')
            os.makedirs(output_file, exist_ok=True)
            output_file_fa = os.path.join(output_file, 'output_mothur.fasta')
            output_file_tax = os.path.join(output_file, 'output_mothur.taxonomy')

            df = read_refcsv(refcsv)
            fasta_mothur = []
            tax_mothur = []
            # 使用name建立联系
            for seq_record in SeqIO.parse(file, 'fasta'):
                name = seq_record.id
                print(name)
                for index, row in df.iterrows():
                    if row['name'] == name:
                        # 处理序列信息
                        kingdom = remove_prefix(row['kingdom'])
                        phylum = remove_prefix(row['phylum'])
                        class_name = remove_prefix(row['class'])
                        order = remove_prefix(row['order'])
                        family = remove_prefix(row['family'])
                        genus = remove_prefix(row['genus'])
                        species = remove_prefix(row['species']).replace(' ', '_')

                        seq_record.id = str(row['name'])
                        seq_record.description = ''
                        fasta_mothur.append(seq_record)

                        new_tax = f"{kingdom};{phylum};{class_name};{order};{family};{genus};{species};"
                        new_tax_df = f"{row['name']}	{new_tax}"
                        tax_mothur.append(new_tax_df)

            SeqIO.write(fasta_mothur, output_file_fa, 'fasta')

            with open(output_file_tax, 'w', encoding='utf_8_sig') as f:
                f.write('\n'.join(tax_mothur))

            print('finish mothur output')

def output(type, file, refcsv):
    if type == 'qiime2':
        output_qiime(fasta_file=file, refcsv=refcsv)
    if type == 'sintax':
        output_sintax(fasta_file=file, refcsv=refcsv)
    if type == 'RDP':
        output_RDP(fasta_file=file, refcsv=refcsv)
    if type == 'BLAST':
        output_BLAST(fasta_file=file, refcsv=refcsv)
    if type == 'kraken2':
        output_kraken2(fasta_file=file, refcsv=refcsv)
    if type == 'DADA2':
        output_DADA2(fasta_file=file, refcsv=refcsv)
    if type == 'mothur':
        output_mothur(fasta_file=file, refcsv=refcsv)

if __name__ == "__main__":
    work_dir = './work'
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    Output_type = 'qiime2'  # qiime2 sintax RDP BLAST kraken2 DADA2 mothur
    output(
        type=Output_type, 
        file=os.path.join(work_dir, 'Filter_5000', 'Filter.fasta'),
        refcsv=os.path.join(work_dir, 'Filter_5000', 'Filter.csv')
    )
