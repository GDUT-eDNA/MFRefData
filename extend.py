###先判断扩展的类型
###再判断扩展的方法（默认/生态学分类等级）
import os
import pandas as pd
from Bio import SeqIO
import numpy as np


prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']


#每个读取函数都会返回fasta和固定格式，用sequenceid建立链接
def read_db2_qiime2(file_path):
    database2_fasta = False
    database2_info_file = False

    for file in os.listdir(file_path):
        file_lower = file.lower()
        if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
            database2_fasta = os.path.join(file_path, file)
        if file_lower.endswith('txt'):
            database2_info_file = os.path.join(file_path, file)

    if database2_fasta == False:
        return print('lacking qiime2 fasta file')
    if database2_info_file == False:
        return print('lacking qiime2 txt file or please ensure taxonomic file end with txt')

    fasta_id = []
    for seq_record in SeqIO.parse(database2_fasta, 'fasta'):
        fasta_id.append(seq_record.id)

    database2_info = [['name','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

    with open(database2_info_file, 'r', encoding='utf_8_sig') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if '\t' not in line:
                continue

            name, tax_str = line.split('\t', 1)

            if name not in fasta_id:
                continue

            items = tax_str.split(';')
            tax_values = ['na', 'na', 'na', 'na', 'na', 'na', 'na']

            for i, pre in enumerate(prefixes):
                for item in items:
                    item = item.strip()
                    if item.startswith(pre):
                        value = item[len(pre):]
                        if value == '':
                            value = 'na'
                        tax_values[i] = value
                        break

            if tax_values[0] != 'na':
                tax_values[0] = 'k_' + tax_values[0]
            if tax_values[1] != 'na':
                tax_values[1] = 'p_' + tax_values[1]
            if tax_values[2] != 'na':
                tax_values[2] = 'c_' + tax_values[2]
            if tax_values[3] != 'na':
                tax_values[3] = 'o_' + tax_values[3]
            if tax_values[4] != 'na':
                tax_values[4] = 'f_' + tax_values[4]
            if tax_values[5] != 'na':
                tax_values[5] = 'g_' + tax_values[5]
            if tax_values[6] != 'na':
                tax_values[6] = 's_' + tax_values[6].replace('_', ' ')

            database2_info.append([name] + tax_values)

    return database2_info, database2_fasta


def read_db2_sintax(file_path):
    database2_fasta = False

    for file in os.listdir(file_path):
        file_lower = file.lower()
        if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
            database2_fasta = os.path.join(file_path, file)

    if database2_fasta == False:
        return print('lacking sintax fasta file')

    database2_info = [['name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

    for seq_record in SeqIO.parse(database2_fasta, 'fasta'):
        header = seq_record.description

        if ';tax=' not in header:
            continue

        name = header.split(';tax=', 1)[0]
        tax_str = header.split(';tax=', 1)[1].strip(';')
        tax_items = tax_str.split(',')

        kingdom = 'na'
        phylum = 'na'
        class_name = 'na'
        order = 'na'
        family = 'na'
        genus = 'na'
        species = 'na'

        for item in tax_items:
            if ':' in item:
                rank, value = item.split(':', 1)
                value = value.strip()

                if rank == 'k' and value != '':
                    if value.startswith('k_'):
                        kingdom = value
                    else:
                        kingdom = 'k_' + value

                if rank == 'p' and value != '':
                    if value.startswith('p_'):
                        phylum = value
                    else:
                        phylum = 'p_' + value

                if rank == 'c' and value != '':
                    if value.startswith('c_'):
                        class_name = value
                    else:
                        class_name = 'c_' + value

                if rank == 'o' and value != '':
                    if value.startswith('o_'):
                        order = value
                    else:
                        order = 'o_' + value

                if rank == 'f' and value != '':
                    if value.startswith('f_'):
                        family = value
                    else:
                        family = 'f_' + value

                if rank == 'g' and value != '':
                    if value.startswith('g_'):
                        genus = value
                    else:
                        genus = 'g_' + value

                if rank == 's' and value != '':
                    value = value.replace('_', ' ')
                    if value.startswith('s_'):
                        species = value
                    else:
                        species = 's_' + value

        database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

    return database2_info, database2_fasta


def read_db2_RDP(file_path):
    database2_fasta = False

    for file in os.listdir(file_path):
        file_lower = file.lower()
        if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
            database2_fasta = os.path.join(file_path, file)

    if database2_fasta == False:
        return print('lacking RDP fasta file')

    database2_info = [['name','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

    for seq_record in SeqIO.parse(database2_fasta,'fasta'):
        header = seq_record.description
        if ' ' not in header:
            continue

        name = header.split(' ', 1)[0]
        tax_str = header.split(' ', 1)[1]
        tax_items = tax_str.strip(';').split(';')

        if len(tax_items) > 0 and tax_items[0] == 'Root':
            tax_items = tax_items[1:]

        kingdom = 'na'
        phylum = 'na'
        class_name = 'na'
        order = 'na'
        family = 'na'
        genus = 'na'
        species = 'na'

        if len(tax_items) > 0 and tax_items[0] != '':
            kingdom = 'k_' + tax_items[0]
        if len(tax_items) > 1 and tax_items[1] != '':
            phylum = 'p_' + tax_items[1]
        if len(tax_items) > 2 and tax_items[2] != '':
            class_name = 'c_' + tax_items[2]
        if len(tax_items) > 3 and tax_items[3] != '':
            order = 'o_' + tax_items[3]
        if len(tax_items) > 4 and tax_items[4] != '':
            family = 'f_' + tax_items[4]
        if len(tax_items) > 5 and tax_items[5] != '':
            genus = 'g_' + tax_items[5]
        if len(tax_items) > 6 and tax_items[6] != '':
            species = 's_' + tax_items[6].replace('_', ' ')

        database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

    return database2_info, database2_fasta


def read_db2_BLAST(file_path):
    database2_fasta = False
    database2_info_file = False

    for file in os.listdir(file_path):
        file_lower = file.lower()
        if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
            database2_fasta = os.path.join(file_path, file)
        if file_lower.endswith('tsv') or file_lower.endswith('csv'):
            database2_info_file = os.path.join(file_path, file)

    if database2_fasta == False:
        return print('lacking BLAST fasta file')

    fasta_id = []
    for seq_record in SeqIO.parse(database2_fasta, 'fasta'):
        fasta_id.append(seq_record.id)

    database2_info = [['name','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

    if database2_info_file != False:
        try:
            df = pd.read_csv(database2_info_file, sep=None, engine='python', encoding='utf_8_sig')
            df.columns = [str(col).strip() for col in df.columns]

            if 'name' in df.columns:
                has_tax = False
                for col in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'scientific name']:
                    if col in df.columns:
                        has_tax = True

                if has_tax == True:
                    for index, row in df.iterrows():
                        name = str(row['name'])

                        if name not in fasta_id:
                            continue

                        kingdom = 'na'
                        phylum = 'na'
                        class_name = 'na'
                        order = 'na'
                        family = 'na'
                        genus = 'na'
                        species = 'na'

                        if 'kingdom' in df.columns and pd.isna(row['kingdom']) == False:
                            value = str(row['kingdom']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('k_'):
                                    kingdom = value
                                else:
                                    kingdom = 'k_' + value

                        if 'phylum' in df.columns and pd.isna(row['phylum']) == False:
                            value = str(row['phylum']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('p_'):
                                    phylum = value
                                else:
                                    phylum = 'p_' + value

                        if 'class' in df.columns and pd.isna(row['class']) == False:
                            value = str(row['class']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('c_'):
                                    class_name = value
                                else:
                                    class_name = 'c_' + value

                        if 'order' in df.columns and pd.isna(row['order']) == False:
                            value = str(row['order']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('o_'):
                                    order = value
                                else:
                                    order = 'o_' + value

                        if 'family' in df.columns and pd.isna(row['family']) == False:
                            value = str(row['family']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('f_'):
                                    family = value
                                else:
                                    family = 'f_' + value

                        if 'genus' in df.columns and pd.isna(row['genus']) == False:
                            value = str(row['genus']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('g_'):
                                    genus = value
                                else:
                                    genus = 'g_' + value

                        if 'species' in df.columns and pd.isna(row['species']) == False:
                            value = str(row['species']).strip().replace('_', ' ')
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('s_'):
                                    species = value
                                else:
                                    species = 's_' + value

                        if species == 'na' and 'scientific name' in df.columns and pd.isna(row['scientific name']) == False:
                            value = str(row['scientific name']).strip().replace('_', ' ')
                            if value != '' and value.lower() != 'nan':
                                species = 's_' + value
                                genus = 'g_' + value.split(' ')[0]

                        database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

                    if len(database2_info) > 1:
                        return database2_info, database2_fasta

        except:
            pass

    for seq_record in SeqIO.parse(database2_fasta,'fasta'):
        name = seq_record.id
        header = seq_record.description

        kingdom = 'na'
        phylum = 'na'
        class_name = 'na'
        order = 'na'
        family = 'na'
        genus = 'na'
        species = 'na'

        if len(header.split(' ', 1)) > 1:
            species_value = header.split(' ', 1)[1].replace('_', ' ')
            if species_value != '':
                species = 's_' + species_value
                genus = 'g_' + species_value.split(' ')[0]

        database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

    return database2_info, database2_fasta


def read_db2_kraken2(file_path):
    database2_fasta = False
    database2_info_file = False

    for file in os.listdir(file_path):
        file_lower = file.lower()
        if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
            database2_fasta = os.path.join(file_path, file)
        if file_lower.endswith('tsv') or file_lower.endswith('csv'):
            database2_info_file = os.path.join(file_path, file)

    if database2_fasta == False:
        return print('lacking kraken2 fasta file')

    fasta_id = []
    for seq_record in SeqIO.parse(database2_fasta, 'fasta'):
        name = seq_record.id
        if '|kraken:taxid|' in name:
            name = name.split('|kraken:taxid|')[0]
        fasta_id.append(name)

    database2_info = [['name','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

    if database2_info_file != False:
        try:
            df = pd.read_csv(database2_info_file, sep=None, engine='python', encoding='utf_8_sig')
            df.columns = [str(col).strip() for col in df.columns]

            if 'name' in df.columns:
                has_tax = False
                for col in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'scientific name']:
                    if col in df.columns:
                        has_tax = True

                if has_tax == True:
                    for index, row in df.iterrows():
                        name = str(row['name'])

                        if name not in fasta_id:
                            continue

                        kingdom = 'na'
                        phylum = 'na'
                        class_name = 'na'
                        order = 'na'
                        family = 'na'
                        genus = 'na'
                        species = 'na'

                        if 'kingdom' in df.columns and pd.isna(row['kingdom']) == False:
                            value = str(row['kingdom']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('k_'):
                                    kingdom = value
                                else:
                                    kingdom = 'k_' + value

                        if 'phylum' in df.columns and pd.isna(row['phylum']) == False:
                            value = str(row['phylum']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('p_'):
                                    phylum = value
                                else:
                                    phylum = 'p_' + value

                        if 'class' in df.columns and pd.isna(row['class']) == False:
                            value = str(row['class']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('c_'):
                                    class_name = value
                                else:
                                    class_name = 'c_' + value

                        if 'order' in df.columns and pd.isna(row['order']) == False:
                            value = str(row['order']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('o_'):
                                    order = value
                                else:
                                    order = 'o_' + value

                        if 'family' in df.columns and pd.isna(row['family']) == False:
                            value = str(row['family']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('f_'):
                                    family = value
                                else:
                                    family = 'f_' + value

                        if 'genus' in df.columns and pd.isna(row['genus']) == False:
                            value = str(row['genus']).strip()
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('g_'):
                                    genus = value
                                else:
                                    genus = 'g_' + value

                        if 'species' in df.columns and pd.isna(row['species']) == False:
                            value = str(row['species']).strip().replace('_', ' ')
                            if value != '' and value.lower() != 'nan':
                                if value.startswith('s_'):
                                    species = value
                                else:
                                    species = 's_' + value

                        if species == 'na' and 'scientific name' in df.columns and pd.isna(row['scientific name']) == False:
                            value = str(row['scientific name']).strip().replace('_', ' ')
                            if value != '' and value.lower() != 'nan':
                                species = 's_' + value
                                genus = 'g_' + value.split(' ')[0]

                        database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

                    if len(database2_info) > 1:
                        return database2_info, database2_fasta

        except:
            pass

    for seq_record in SeqIO.parse(database2_fasta,'fasta'):
        header = seq_record.description
        name = seq_record.id

        if '|kraken:taxid|' in name:
            name = name.split('|kraken:taxid|')[0]

        kingdom = 'na'
        phylum = 'na'
        class_name = 'na'
        order = 'na'
        family = 'na'
        genus = 'na'
        species = 'na'

        if len(header.split(' ', 1)) > 1:
            species_value = header.split(' ', 1)[1].replace('_', ' ')
            if species_value != '':
                species = 's_' + species_value
                genus = 'g_' + species_value.split(' ')[0]

        database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

    return database2_info, database2_fasta




def read_db2_mothur(file_path):
    database2_fasta = False
    database2_info_file = False

    for file in os.listdir(file_path):
        file_lower = file.lower()
        if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
            database2_fasta = os.path.join(file_path, file)
        if file_lower.endswith('taxonomy') or file_lower.endswith('tax'):
            database2_info_file = os.path.join(file_path, file)

    if database2_fasta == False:
        return print('lacking mothur fasta file')
    if database2_info_file == False:
        return print('lacking mothur taxonomy file')

    fasta_id = []
    for seq_record in SeqIO.parse(database2_fasta, 'fasta'):
        fasta_id.append(seq_record.id)

    database2_info = [['name','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

    with open(database2_info_file, 'r', encoding='utf_8_sig') as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if '\t' not in line:
                continue

            name, tax_str = line.split('\t', 1)

            if name not in fasta_id:
                continue

            tax_items = tax_str.strip(';').split(';')

            kingdom = 'na'
            phylum = 'na'
            class_name = 'na'
            order = 'na'
            family = 'na'
            genus = 'na'
            species = 'na'

            if len(tax_items) > 0 and tax_items[0] != '':
                kingdom = 'k_' + tax_items[0]
            if len(tax_items) > 1 and tax_items[1] != '':
                phylum = 'p_' + tax_items[1]
            if len(tax_items) > 2 and tax_items[2] != '':
                class_name = 'c_' + tax_items[2]
            if len(tax_items) > 3 and tax_items[3] != '':
                order = 'o_' + tax_items[3]
            if len(tax_items) > 4 and tax_items[4] != '':
                family = 'f_' + tax_items[4]
            if len(tax_items) > 5 and tax_items[5] != '':
                genus = 'g_' + tax_items[5]
            if len(tax_items) > 6 and tax_items[6] != '':
                species = 's_' + tax_items[6].replace('_', ' ')

            database2_info.append([name, kingdom, phylum, class_name, order, family, genus, species])

    return database2_info, database2_fasta


def extend(formation,db1,db2,level = False,value = False,output = False):
    print('module Extend for combining two reference database')
    print('please ensure MFrefdata database file contain csv file and combined fasta file')
    print('AND make sure db1 type is MF')
    print('\n')
    print('the other reference database type require please read')
    print('MFRefData recommend fasta with sequence.id')
    print('\n')
    print('MF means the formation of two reference databases both are MF')

    if db1 and db2:
        if os.path.exists(db1) and os.path.exists(db2):
            #读取db1
            for file in os.listdir(db1):
                file_lower = file.lower()
                if file_lower.endswith('fasta') or file_lower.endswith('fa') or file_lower.endswith('fna'):
                    db1_fasta = os.path.join(db1, file)
                if file_lower.endswith('csv'):
                    db1_df = os.path.join(db1, file)

            df1 = pd.read_csv(db1_df, encoding='utf_8_sig')

            for col in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                if col not in df1.columns:
                    df1[col] = 'na'
                df1[col] = df1[col].fillna('na')

            #根据db2按照格式读取
            if formation == 'qiime2':
                db2_df, db2_fasta = read_db2_qiime2(db2)
            if formation == 'sintax':
                db2_df, db2_fasta = read_db2_sintax(db2)
            if formation == 'RDP':
                db2_df, db2_fasta = read_db2_RDP(db2)
            if formation == 'BLAST':
                db2_df, db2_fasta = read_db2_BLAST(db2)
            if formation == 'kraken2':
                db2_df, db2_fasta = read_db2_kraken2(db2)
            if formation == 'mothur':
                db2_df, db2_fasta = read_db2_mothur(db2)

            db2_df_head = db2_df[0]
            db2_df_info = db2_df[1:]
            df2 = pd.DataFrame(db2_df_info,columns=db2_df_head)

            for col in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                if col not in df2.columns:
                    df2[col] = 'na'
                df2[col] = df2[col].fillna('na')

        else:
            return print('error')

        #输出目录
        if output == False:
            output = os.path.commonpath([db1, db2])

        fp = os.path.join(output,f'Extend_{formation}')
        os.makedirs(fp, exist_ok=True)

        #默认补充
        if level == False:
            name_toremove = []
            new_database = []

            names1 = set(df1['name'])
            names2 = set(df2['name'])
            duplicate_names = names1 & names2
            name_toremove = set(duplicate_names)

            if name_toremove:
                df1 = df1[~df1['name'].isin(name_toremove)]

                for seq_record in SeqIO.parse(db1_fasta,'fasta'):
                    if seq_record.id in name_toremove:
                        continue
                    else:
                        new_database.append(seq_record)
            else:
                for seq_record in SeqIO.parse(db1_fasta,'fasta'):
                    new_database.append(seq_record)

            df2_names = set(df2['name'])

            for seq_record in SeqIO.parse(db2_fasta,'fasta'):
                name = seq_record.id

                if formation == 'sintax':
                    if ';tax=' in name:
                        name = name.split(';tax=')[0]

                if formation == 'kraken2':
                    if '|kraken:taxid|' in name:
                        name = name.split('|kraken:taxid|')[0]

                if name in df2_names:
                    seq_record.id = name
                    seq_record.name = name
                    seq_record.description = ''
                    new_database.append(seq_record)

            SeqIO.write(new_database,os.path.join(fp,'Extend.fasta'),'fasta')

            new_df = pd.concat([df1, df2], ignore_index=True)
            new_df.to_csv(os.path.join(fp,'Extend.csv'), encoding='utf-8', index=False)

            print('finish Extend')
            return

        if level != False:
            value_list = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

            if level in value_list:
                print(value_list)
                print(f'{level} level mapping')

                mapping_name = []

                df2_value = df2[f'{level}'].value_counts().index.tolist()
                set1 = set(df1['name'])
                set2 = set(df2['name'])
                set_to_add = set1 - set2

                if not set_to_add:
                    print('same db, no results')
                    return
                else:
                    #首先判断是否同分类学等级，其次判断是否重复
                    for name in set_to_add:
                        temp = df1.loc[df1['name'] == name,f'{level}'].iloc[0]

                        if temp in df2_value:
                            mapping_name.append(name)

                    df1_add = df1[df1['name'].isin(mapping_name)]

                    new_df = pd.concat([df1_add, df2], ignore_index=True)
                    new_df.to_csv(os.path.join(fp,'Extend.csv'), encoding='utf-8', index=False)

                    new_database = []

                    df2_names = set(df2['name'])

                    for seq_record in SeqIO.parse(db2_fasta,'fasta'):
                        name = seq_record.id

                        if formation == 'sintax':
                            if ';tax=' in name:
                                name = name.split(';tax=')[0]

                        if formation == 'kraken2':
                            if '|kraken:taxid|' in name:
                                name = name.split('|kraken:taxid|')[0]

                        if name in df2_names:
                            seq_record.id = name
                            seq_record.name = name
                            seq_record.description = ''
                            new_database.append(seq_record)

                    for seq_record in SeqIO.parse(db1_fasta,'fasta'):
                        if seq_record.id in mapping_name:
                            new_database.append(seq_record)

                    SeqIO.write(new_database, os.path.join(fp,'Extend.fasta'), 'fasta')

                    print('finish Extend')
                    print(f'add {len(mapping_name)} sequences to {db2}')
                    return

            else:
                return print('please check level')

    return print('finish extend')
