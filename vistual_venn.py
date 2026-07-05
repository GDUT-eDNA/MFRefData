import os
import re
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from matplotlib.patches import Patch
from matplotlib_venn import venn2, venn3
from matplotlib_venn.layout.venn2 import DefaultLayoutAlgorithm


# =========================
# 1. 全局参数
# =========================

RANK_COLUMNS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

FASTA_SUFFIX = ('fasta', 'fa', 'fna')
TABLE_SUFFIX = ('csv', 'tsv')
QIIME2_TAX_SUFFIX = ('txt',)
MOTHUR_TAX_SUFFIX = ('taxonomy', 'tax')

PLOT_PARAMS = {
    'dpi': 400,

    'venn2_figsize': (8, 5),
    'venn2_normalize_to': 2.0,
    'venn2_xlim': (-2.0, 2.0),
    'venn2_ylim': (-1.8, 1.8),
    'venn2_colors': ('steelblue', 'gold'),
    'venn2_alpha': 0.7,

    'venn3_figsize': (24, 16),
    'venn3_colors': ('steelblue', 'red', 'green'),
    'venn3_alpha': 0.5,

    'title_fontsize': 15,
    'legend_fontsize': 8,
    'text_fontsize': 7,
    'common_text_fontsize': 6,

    'edge_color': 'black',
    'edge_width': 2,

    'venn3_text_positions': {
        'only_species_list': (-0.95, 0.5),
        'only_mf': (-0.95, 0),
        'only_third': (1.1, 0.5),
        'common': (-0.95, -0.62),
        'species_third_common': (1.1, -0.12)
    },

    'venn3_text_colors': {
        'only_species_list': 'steelblue',
        'only_mf': 'red',
        'only_third': 'green',
        'common': 'purple',
        'species_third_common': 'gold'
    }
}


# =========================
# 2. 文件自动检索
# =========================

def find_files_by_suffix(file_path, suffix_tuple):
    result_files = []

    if os.path.isfile(file_path):
        file_lower = os.path.basename(file_path).lower()
        if file_lower.endswith(suffix_tuple):
            result_files.append(file_path)

    if os.path.isdir(file_path):
        for file in os.listdir(file_path):
            file_lower = file.lower()
            if file_lower.endswith(suffix_tuple):
                result_files.append(os.path.join(file_path, file))

    return result_files


def find_one_file(file_path, suffix_tuple):
    result_files = find_files_by_suffix(file_path, suffix_tuple)

    if len(result_files) == 0:
        return False

    return result_files[0]


def find_file_by_keyword(file_path, suffix_tuple, keyword):
    result_files = find_files_by_suffix(file_path, suffix_tuple)

    for file in result_files:
        file_lower = os.path.basename(file).lower()
        if keyword.lower() in file_lower:
            return file

    return False


def check_required_files(formation, format_path):
    formation = str(formation).upper()

    if os.path.exists(format_path) == False:
        print('format path is not exist')
        return False

    if formation == 'QIIME2':
        fasta_file = find_one_file(format_path, FASTA_SUFFIX)
        txt_file = find_one_file(format_path, QIIME2_TAX_SUFFIX)

        if fasta_file == False:
            print('lacking QIIME2 fasta file')
            return False
        if txt_file == False:
            print('lacking QIIME2 taxonomy txt file')
            return False

        return {'fasta': fasta_file, 'txt': txt_file}

    if formation == 'SINTAX':
        fasta_file = find_one_file(format_path, FASTA_SUFFIX)

        if fasta_file == False:
            print('lacking SINTAX fasta file')
            return False

        return {'fasta': fasta_file}

    if formation == 'RDP':
        fasta_file = find_one_file(format_path, FASTA_SUFFIX)

        if fasta_file == False:
            print('lacking RDP fasta file')
            return False

        return {'fasta': fasta_file}

    if formation == 'BLAST':
        fasta_file = find_one_file(format_path, FASTA_SUFFIX)
        table_file = find_one_file(format_path, TABLE_SUFFIX)

        if fasta_file == False:
            print('lacking BLAST fasta file')
            return False
        if table_file == False:
            print('lacking BLAST tsv/csv file')
            return False

        return {'fasta': fasta_file, 'table': table_file}

    if formation == 'KRAKEN2':
        fasta_file = find_one_file(format_path, FASTA_SUFFIX)

        if fasta_file == False:
            print('lacking Kraken2 fasta file')
            return False

        return {'fasta': fasta_file}

    if formation == 'DADA2':
        full_file = find_file_by_keyword(format_path, FASTA_SUFFIX, 'full')

        if full_file == False:
            full_file = find_file_by_keyword(format_path, FASTA_SUFFIX, 'taxonomy')

        species_file = find_file_by_keyword(format_path, FASTA_SUFFIX, 'genus_species')

        if species_file == False:
            species_file = find_file_by_keyword(format_path, FASTA_SUFFIX, 'species')

        if full_file == False and species_file == False:
            print('lacking DADA2 fasta file')
            return False

        return {'full': full_file, 'species': species_file}

    if formation == 'MOTHUR':
        fasta_file = find_one_file(format_path, FASTA_SUFFIX)
        taxonomy_file = find_one_file(format_path, MOTHUR_TAX_SUFFIX)

        if fasta_file == False:
            print('lacking mothur fasta file')
            return False
        if taxonomy_file == False:
            print('lacking mothur taxonomy file')
            return False

        return {'fasta': fasta_file, 'taxonomy': taxonomy_file}

    print('unsupported formation')
    return False


# =========================
# 3. 基础处理函数
# =========================

def empty_tax_row():
    return {
        'kingdom': 'na',
        'phylum': 'na',
        'class': 'na',
        'order': 'na',
        'family': 'na',
        'genus': 'na',
        'species': 'na'
    }


def normalize_empty(value):
    if pd.isna(value):
        return 'na'

    value = str(value).strip()

    if value == '' or value.lower() == 'nan' or value.lower() == 'na':
        return 'na'

    return value


def normalize_species(value):
    value = normalize_empty(value)

    if value == 'na':
        return 'na'

    value = value.replace('_', ' ')
    value = re.sub(r'\s+', ' ', value).strip()

    if value == '':
        return 'na'

    return value


def build_taxonomy_df(result_list):
    taxonomy_df = pd.DataFrame(result_list, columns=RANK_COLUMNS)
    taxonomy_df = taxonomy_df.drop_duplicates()
    return taxonomy_df


def build_unique_rank_table(taxonomy_df):
    unique_rank_dict = {}
    max_len = 0

    for rank in RANK_COLUMNS:
        values = []

        for value in taxonomy_df[rank].tolist():
            if rank == 'species':
                value = normalize_species(value)
            else:
                value = normalize_empty(value)

            if value != 'na' and value not in values:
                values.append(value)

        unique_rank_dict[rank] = values

        if len(values) > max_len:
            max_len = len(values)

    for rank in RANK_COLUMNS:
        while len(unique_rank_dict[rank]) < max_len:
            unique_rank_dict[rank].append('')

    unique_rank_df = pd.DataFrame(unique_rank_dict, columns=RANK_COLUMNS)

    return unique_rank_df


def save_species_detail_table(group_dict, output_file):
    result_list = []

    for group_name, species_set in group_dict.items():
        for species in sorted(list(species_set)):
            result_list.append({
                'group': group_name,
                'species': species
            })

    detail_df = pd.DataFrame(result_list, columns=['group', 'species'])
    detail_df.to_csv(output_file, index=False, encoding='utf-8-sig')

    print(f'common/different species saved to: {output_file}')

    return detail_df


# =========================
# 4. 读取物种名录与 MF 数据库
# =========================

def read_species_list(species_list):
    print('reading species list')

    if os.path.exists(species_list) == False:
        print('species list path is not exist')
        return False

    try:
        list_df = pd.read_csv(species_list, sep='\t', header=None, encoding='utf-8-sig')
    except Exception as e:
        print(f'Error reading species list file: {e}')
        return False

    species_set = set()

    for value in list_df[0].value_counts().index.tolist():
        species = normalize_species(value)

        if species != 'na':
            species_set.add(species)

    return species_set


def read_MF_database(file):
    print('reading MF database')

    if os.path.exists(file) == False:
        print('MF csv path is not exist')
        return False

    try:
        csv_df = pd.read_csv(file, encoding='utf-8-sig')
    except Exception as e:
        print(f'Error reading MF csv file: {e}')
        return False

    result_list = []

    for index, row in csv_df.iterrows():
        result = empty_tax_row()

        if 'kingdom' in csv_df.columns:
            value = normalize_empty(row['kingdom'])
            if value.startswith('k_'):
                value = value[2:]
            result['kingdom'] = normalize_empty(value)

        if 'phylum' in csv_df.columns:
            value = normalize_empty(row['phylum'])
            if value.startswith('p_'):
                value = value[2:]
            result['phylum'] = normalize_empty(value)

        if 'class' in csv_df.columns:
            value = normalize_empty(row['class'])
            if value.startswith('c_'):
                value = value[2:]
            result['class'] = normalize_empty(value)

        if 'order' in csv_df.columns:
            value = normalize_empty(row['order'])
            if value.startswith('o_'):
                value = value[2:]
            result['order'] = normalize_empty(value)

        if 'family' in csv_df.columns:
            value = normalize_empty(row['family'])
            if value.startswith('f_'):
                value = value[2:]
            result['family'] = normalize_empty(value)

        if 'genus' in csv_df.columns:
            value = normalize_empty(row['genus'])
            if value.startswith('g_'):
                value = value[2:]
            result['genus'] = normalize_empty(value)

        if 'species' in csv_df.columns:
            value = normalize_empty(row['species'])
            if value.startswith('s_'):
                value = value[2:]
            result['species'] = normalize_species(value)

        if result['species'] == 'na' and 'scientific name' in csv_df.columns:
            result['species'] = normalize_species(row['scientific name'])

        result_list.append(result)

    return build_taxonomy_df(result_list)


# =========================
# 5. 分格式读取第三方数据库
# =========================

def read_qiime2_database(format_path):
    print('reading QIIME2 database')

    files = check_required_files('QIIME2', format_path)

    if files == False:
        return False

    fasta_id = []
    for seq_record in SeqIO.parse(files['fasta'], 'fasta'):
        fasta_id.append(seq_record.id)

    result_list = []

    with open(files['txt'], 'r', encoding='utf-8-sig') as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if '\t' not in line:
                continue

            name, tax_str = line.split('\t', 1)

            if name not in fasta_id:
                continue

            tax_items = tax_str.split(';')
            result = empty_tax_row()

            for item in tax_items:
                item = item.strip()

                if item.startswith('k__'):
                    result['kingdom'] = normalize_empty(item[3:])

                if item.startswith('p__'):
                    result['phylum'] = normalize_empty(item[3:])

                if item.startswith('c__'):
                    result['class'] = normalize_empty(item[3:])

                if item.startswith('o__'):
                    result['order'] = normalize_empty(item[3:])

                if item.startswith('f__'):
                    result['family'] = normalize_empty(item[3:])

                if item.startswith('g__'):
                    result['genus'] = normalize_empty(item[3:])

                if item.startswith('s__'):
                    result['species'] = normalize_species(item[3:])

            result_list.append(result)

    return build_taxonomy_df(result_list)


def read_sintax_database(format_path):
    print('reading SINTAX database')

    files = check_required_files('SINTAX', format_path)

    if files == False:
        return False

    result_list = []

    for seq_record in SeqIO.parse(files['fasta'], 'fasta'):
        header = seq_record.description

        if ';tax=' not in header:
            continue

        tax_str = header.split(';tax=', 1)[1].strip(';')
        tax_items = tax_str.split(',')
        result = empty_tax_row()

        for item in tax_items:
            item = item.strip()

            if ':' not in item:
                continue

            rank, value = item.split(':', 1)
            value = normalize_empty(value)

            if rank == 'k':
                if value.startswith('k_'):
                    value = value[2:]
                result['kingdom'] = normalize_empty(value)

            if rank == 'p':
                if value.startswith('p_'):
                    value = value[2:]
                result['phylum'] = normalize_empty(value)

            if rank == 'c':
                if value.startswith('c_'):
                    value = value[2:]
                result['class'] = normalize_empty(value)

            if rank == 'o':
                if value.startswith('o_'):
                    value = value[2:]
                result['order'] = normalize_empty(value)

            if rank == 'f':
                if value.startswith('f_'):
                    value = value[2:]
                result['family'] = normalize_empty(value)

            if rank == 'g':
                if value.startswith('g_'):
                    value = value[2:]
                result['genus'] = normalize_empty(value)

            if rank == 's':
                if value.startswith('s_'):
                    value = value[2:]
                result['species'] = normalize_species(value)

        result_list.append(result)

    return build_taxonomy_df(result_list)


def read_RDP_database(format_path):
    print('reading RDP database')

    files = check_required_files('RDP', format_path)

    if files == False:
        return False

    result_list = []

    for seq_record in SeqIO.parse(files['fasta'], 'fasta'):
        header = seq_record.description

        if ' ' not in header:
            continue

        tax_str = header.split(' ', 1)[1]
        tax_items = tax_str.strip(';').split(';')

        if len(tax_items) > 0 and tax_items[0] == 'Root':
            tax_items = tax_items[1:]

        result = empty_tax_row()

        if len(tax_items) > 0:
            result['kingdom'] = normalize_empty(tax_items[0])
        if len(tax_items) > 1:
            result['phylum'] = normalize_empty(tax_items[1])
        if len(tax_items) > 2:
            result['class'] = normalize_empty(tax_items[2])
        if len(tax_items) > 3:
            result['order'] = normalize_empty(tax_items[3])
        if len(tax_items) > 4:
            result['family'] = normalize_empty(tax_items[4])
        if len(tax_items) > 5:
            result['genus'] = normalize_empty(tax_items[5])
        if len(tax_items) > 6:
            result['species'] = normalize_species(tax_items[6])

        result_list.append(result)

    return build_taxonomy_df(result_list)


def read_BLAST_database(format_path):
    print('reading BLAST database')

    files = check_required_files('BLAST', format_path)

    if files == False:
        return False

    result_list = []

    try:
        table_df = pd.read_csv(files['table'], sep=None, engine='python', encoding='utf-8-sig')
    except Exception:
        try:
            table_df = pd.read_csv(files['table'], sep=None, engine='python', header=None, encoding='utf-8-sig')
        except Exception:
            table_df = False

    if table_df is not False:
        if 'name' not in table_df.columns:
            if table_df.shape[1] >= 2:
                table_df = table_df.iloc[:, 0:2]
                table_df.columns = ['name', 'taxid']

        if 'name' in table_df.columns:
            has_taxonomy = False

            for col in RANK_COLUMNS:
                if col in table_df.columns:
                    has_taxonomy = True

            if 'scientific name' in table_df.columns:
                has_taxonomy = True

            if has_taxonomy == True:
                for index, row in table_df.iterrows():
                    result = empty_tax_row()

                    if 'kingdom' in table_df.columns:
                        value = normalize_empty(row['kingdom'])
                        if value.startswith('k_'):
                            value = value[2:]
                        result['kingdom'] = normalize_empty(value)

                    if 'phylum' in table_df.columns:
                        value = normalize_empty(row['phylum'])
                        if value.startswith('p_'):
                            value = value[2:]
                        result['phylum'] = normalize_empty(value)

                    if 'class' in table_df.columns:
                        value = normalize_empty(row['class'])
                        if value.startswith('c_'):
                            value = value[2:]
                        result['class'] = normalize_empty(value)

                    if 'order' in table_df.columns:
                        value = normalize_empty(row['order'])
                        if value.startswith('o_'):
                            value = value[2:]
                        result['order'] = normalize_empty(value)

                    if 'family' in table_df.columns:
                        value = normalize_empty(row['family'])
                        if value.startswith('f_'):
                            value = value[2:]
                        result['family'] = normalize_empty(value)

                    if 'genus' in table_df.columns:
                        value = normalize_empty(row['genus'])
                        if value.startswith('g_'):
                            value = value[2:]
                        result['genus'] = normalize_empty(value)

                    if 'species' in table_df.columns:
                        value = normalize_empty(row['species'])
                        if value.startswith('s_'):
                            value = value[2:]
                        result['species'] = normalize_species(value)

                    if result['species'] == 'na' and 'scientific name' in table_df.columns:
                        result['species'] = normalize_species(row['scientific name'])
                        if result['species'] != 'na':
                            result['genus'] = result['species'].split(' ')[0]

                    result_list.append(result)

                if len(result_list) > 0:
                    return build_taxonomy_df(result_list)

    for seq_record in SeqIO.parse(files['fasta'], 'fasta'):
        result = empty_tax_row()
        header = seq_record.description

        if ' ' in header:
            species_info = header.split(' ', 1)[1]
            species_info = normalize_species(species_info)

            species_items = species_info.split()

            if len(species_items) >= 2:
                result['genus'] = species_items[0]
                result['species'] = species_items[0] + ' ' + species_items[1]

        result_list.append(result)

    return build_taxonomy_df(result_list)


def read_kraken2_database(format_path):
    print('reading Kraken2 database')

    files = check_required_files('KRAKEN2', format_path)

    if files == False:
        return False

    result_list = []

    for seq_record in SeqIO.parse(files['fasta'], 'fasta'):
        result = empty_tax_row()
        header = seq_record.description

        if ' ' in header:
            species_info = header.split(' ', 1)[1]
            species_info = normalize_species(species_info)

            species_items = species_info.split()

            if len(species_items) >= 2:
                result['genus'] = species_items[0]
                result['species'] = species_items[0] + ' ' + species_items[1]

        result_list.append(result)

    return build_taxonomy_df(result_list)


def read_DADA2_database(format_path):
    print('reading DADA2 database')

    files = check_required_files('DADA2', format_path)

    if files == False:
        return False

    result_list = []

    if files['full'] != False:
        for seq_record in SeqIO.parse(files['full'], 'fasta'):
            tax_items = seq_record.description.strip(';').split(';')
            result = empty_tax_row()

            if len(tax_items) > 0:
                result['kingdom'] = normalize_empty(tax_items[0])
            if len(tax_items) > 1:
                result['phylum'] = normalize_empty(tax_items[1])
            if len(tax_items) > 2:
                result['class'] = normalize_empty(tax_items[2])
            if len(tax_items) > 3:
                result['order'] = normalize_empty(tax_items[3])
            if len(tax_items) > 4:
                result['family'] = normalize_empty(tax_items[4])
            if len(tax_items) > 5:
                result['genus'] = normalize_empty(tax_items[5])
            if len(tax_items) > 6:
                result['species'] = normalize_species(tax_items[6])

            result_list.append(result)

    if len(result_list) == 0 and files['species'] != False:
        for seq_record in SeqIO.parse(files['species'], 'fasta'):
            header_items = seq_record.description.split()
            result = empty_tax_row()

            if len(header_items) >= 3:
                result['genus'] = normalize_empty(header_items[1])
                result['species'] = normalize_species(header_items[1] + ' ' + header_items[2])

            result_list.append(result)

    if len(result_list) == 0:
        print('DADA2 format file cannot be parsed')
        return False

    return build_taxonomy_df(result_list)


def read_mothur_database(format_path):
    print('reading mothur database')

    files = check_required_files('MOTHUR', format_path)

    if files == False:
        return False

    fasta_id = []

    for seq_record in SeqIO.parse(files['fasta'], 'fasta'):
        fasta_id.append(seq_record.id)

    result_list = []

    with open(files['taxonomy'], 'r', encoding='utf-8-sig') as f:
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
            result = empty_tax_row()

            if len(tax_items) > 0:
                value = normalize_empty(tax_items[0])
                if value.startswith('k:'):
                    value = value[2:]
                result['kingdom'] = normalize_empty(value)

            if len(tax_items) > 1:
                value = normalize_empty(tax_items[1])
                if value.startswith('p:'):
                    value = value[2:]
                result['phylum'] = normalize_empty(value)

            if len(tax_items) > 2:
                value = normalize_empty(tax_items[2])
                if value.startswith('c:'):
                    value = value[2:]
                result['class'] = normalize_empty(value)

            if len(tax_items) > 3:
                value = normalize_empty(tax_items[3])
                if value.startswith('o:'):
                    value = value[2:]
                result['order'] = normalize_empty(value)

            if len(tax_items) > 4:
                value = normalize_empty(tax_items[4])
                if value.startswith('f:'):
                    value = value[2:]
                result['family'] = normalize_empty(value)

            if len(tax_items) > 5:
                value = normalize_empty(tax_items[5])
                if value.startswith('g:'):
                    value = value[2:]
                result['genus'] = normalize_empty(value)

            if len(tax_items) > 6:
                value = normalize_empty(tax_items[6])
                if value.startswith('s:'):
                    value = value[2:]
                result['species'] = normalize_species(value)

            result_list.append(result)

    return build_taxonomy_df(result_list)


def read_third_database(formation, format_path):
    formation = str(formation).upper()

    if formation == 'QIIME2':
        return read_qiime2_database(format_path)

    if formation == 'SINTAX':
        return read_sintax_database(format_path)

    if formation == 'RDP':
        return read_RDP_database(format_path)

    if formation == 'BLAST':
        return read_BLAST_database(format_path)

    if formation == 'KRAKEN2':
        return read_kraken2_database(format_path)

    if formation == 'DADA2':
        return read_DADA2_database(format_path)

    if formation == 'MOTHUR':
        return read_mothur_database(format_path)

    print('unsupported formation')
    return False


# =========================
# 6. 集合统计
# =========================

def get_rank_set(tax_df, rank):
    rank_set = set()

    for value in tax_df[rank].value_counts().index.tolist():
        if rank == 'species':
            value = normalize_species(value)
        else:
            value = normalize_empty(value)

        if value != 'na':
            rank_set.add(value)

    return rank_set


def compare_taxonomy_overlap(mf_df, third_df):
    result_list = []

    for rank in RANK_COLUMNS:
        mf_set = get_rank_set(mf_df, rank)
        third_set = get_rank_set(third_df, rank)

        common_set = mf_set & third_set
        only_mf = mf_set - third_set
        only_third = third_set - mf_set

        result_list.append({
            'rank': rank,
            'MF_count': len(mf_set),
            'third_count': len(third_set),
            'common_count': len(common_set),
            'only_MF_count': len(only_mf),
            'only_third_count': len(only_third),
            'common_taxa': '; '.join(sorted(list(common_set))),
            'only_MF_taxa': '; '.join(sorted(list(only_mf))),
            'only_third_taxa': '; '.join(sorted(list(only_third)))
        })

    return pd.DataFrame(result_list)


def make_text(title, value_set, n=2):
    value_list = sorted(list(value_set))[:30]
    text = [title + '\n']

    for i in range(0, len(value_list), n):
        text.extend(value_list[i:i + n])
        text.append('\n')

    return '  '.join(text)


# =========================
# 7. 绘图函数
# =========================

def draw_venn2(species_list_set, mf_species_set, output_dir):
    common = species_list_set & mf_species_set
    only_1 = species_list_set - mf_species_set
    only_2 = mf_species_set - species_list_set

    plt.figure(figsize=PLOT_PARAMS['venn2_figsize'])

    plt_venn2 = venn2(
        [species_list_set, mf_species_set],
        layout_algorithm=DefaultLayoutAlgorithm(normalize_to=PLOT_PARAMS['venn2_normalize_to']),
        set_labels=(
            f'species list\n{len(species_list_set)} species',
            f'MF reference DB\n{len(mf_species_set)} species'
        ),
        set_colors=PLOT_PARAMS['venn2_colors'],
        alpha=PLOT_PARAMS['venn2_alpha'],
        subset_label_formatter=lambda x: ''
    )

    ax = plt.gca()
    ax.set_xlim(*PLOT_PARAMS['venn2_xlim'])
    ax.set_ylim(*PLOT_PARAMS['venn2_ylim'])
    ax.set_title(
        'Species Composition in Species List and MF Reference Database',
        fontsize=PLOT_PARAMS['title_fontsize'],
        y=1.02
    )

    for patch in plt_venn2.patches:
        if patch is not None:
            patch.set_edgecolor(PLOT_PARAMS['edge_color'])
            patch.set_linewidth(PLOT_PARAMS['edge_width'])
            patch.set_linestyle('-')

    color_10 = plt_venn2.get_patch_by_id('10').get_facecolor() if plt_venn2.get_patch_by_id('10') else 'steelblue'
    color_01 = plt_venn2.get_patch_by_id('01').get_facecolor() if plt_venn2.get_patch_by_id('01') else 'gold'
    color_11 = plt_venn2.get_patch_by_id('11').get_facecolor() if plt_venn2.get_patch_by_id('11') else 'gray'

    legend_elements = [
        Patch(facecolor=color_10, label=f'Only in Species List: {len(only_1)}'),
        Patch(facecolor=color_01, label=f'Only in MF Reference DB: {len(only_2)}'),
        Patch(facecolor=color_11, label=f'Both in Two Datasets: {len(common)}')
    ]

    ax.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(0.5, 0.95),
        fontsize=PLOT_PARAMS['legend_fontsize'],
        frameon=True
    )

    common_text = make_text('species in both datasets (top 30):', common, n=5)

    ax.text(
        0.01, 0.12, common_text,
        transform=ax.transAxes,
        fontsize=PLOT_PARAMS['common_text_fontsize'],
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor=color_11, alpha=0.4)
    )

    output_png = os.path.join(output_dir, 'compare_venn2_species.png')

    plt.tight_layout()
    plt.savefig(output_png, format='png', dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
    plt.close()

    print(f'finish venn2: {output_png}')


def draw_venn3(species_list_set, mf_species_set, third_species_set, formation, output_dir):
    common123 = species_list_set & mf_species_set & third_species_set
    both13 = species_list_set & third_species_set

    only_1 = species_list_set - mf_species_set - third_species_set
    only_2 = mf_species_set - species_list_set - third_species_set
    only_3 = third_species_set - species_list_set - mf_species_set

    plt.figure(figsize=PLOT_PARAMS['venn3_figsize'])
    plt.tight_layout(rect=[0.03, 0.3, 0.8, 0.8])

    plt_venn3 = venn3(
        [species_list_set, mf_species_set, third_species_set],
        set_labels=['species list', 'MF reference DB', f'{formation} DB'],
        set_colors=PLOT_PARAMS['venn3_colors'],
        alpha=PLOT_PARAMS['venn3_alpha']
    )

    plt.title('Species Composition in Species List, MF Database and Third-party Database')

    legend_elements = [
        Patch(facecolor='steelblue', label=f'Only in Species List: {len(only_1)}'),
        Patch(facecolor='red', label=f'Only in MF Reference DB: {len(only_2)}'),
        Patch(facecolor='green', label=f'Only in {formation} DB: {len(only_3)}'),
        Patch(facecolor='purple', label=f'Common species in three datasets: {len(common123)}'),
        Patch(facecolor='gold', label=f'Both in Species List and {formation} DB: {len(both13)}')
    ]

    plt.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.66, 0.1),
        fontsize=PLOT_PARAMS['legend_fontsize'],
        ncol=3,
        frameon=True
    )

    legend_text1 = make_text('species only in species list:', only_1, n=3)
    legend_text2 = make_text('species only in MF reference DB:', only_2, n=3)
    legend_text3 = make_text(f'species only in {formation} DB:', only_3, n=3)
    legend_text4 = make_text('species in common:', common123, n=10)
    legend_text5 = make_text(f'Both in species list and {formation} DB:', both13, n=3)

    pos = PLOT_PARAMS['venn3_text_positions']
    colors = PLOT_PARAMS['venn3_text_colors']

    plt.text(*pos['only_species_list'], legend_text1, fontsize=PLOT_PARAMS['text_fontsize'], verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['only_species_list'], alpha=0.4))
    plt.text(*pos['only_mf'], legend_text2, fontsize=PLOT_PARAMS['text_fontsize'], verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['only_mf'], alpha=0.4))
    plt.text(*pos['only_third'], legend_text3, fontsize=PLOT_PARAMS['text_fontsize'], verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['only_third'], alpha=0.4))
    plt.text(*pos['common'], legend_text4, fontsize=PLOT_PARAMS['text_fontsize'], verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['common'], alpha=0.4))
    plt.text(*pos['species_third_common'], legend_text5, fontsize=PLOT_PARAMS['text_fontsize'], verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['species_third_common'], alpha=0.4))

    venn_regions = ['100', '010', '001', '110', '101', '011', '111']

    for region_id in venn_regions:
        patch = plt_venn3.get_patch_by_id(region_id)
        if patch is not None:
            patch.set_edgecolor(PLOT_PARAMS['edge_color'])
            patch.set_linewidth(PLOT_PARAMS['edge_width'])
            patch.set_linestyle('-')

    output_png = os.path.join(output_dir, f'compare_{formation}_venn3_species.png')

    plt.savefig(output_png, format='png', dpi=PLOT_PARAMS['dpi'], bbox_inches='tight')
    plt.close()

    print(f'finish venn3: {output_png}')


# =========================
# 8. 主函数
# =========================

def venn(file, species_list, formation=False, format_path=False, output_dir=False):
    if output_dir == False:
        output_dir = os.path.dirname(file)

    os.makedirs(output_dir, exist_ok=True)

    species_list_set = read_species_list(species_list)

    if species_list_set is False:
        return

    mf_df = read_MF_database(file)

    if mf_df is False:
        return

    mf_unique_df = build_unique_rank_table(mf_df)

    mf_unique_df.to_csv(
        os.path.join(output_dir, 'MF_standard_taxonomy_table.csv'),
        index=False,
        encoding='utf-8-sig'
    )

    mf_species_set = get_rank_set(mf_df, 'species')

    if formation == False:
        species_detail_file = os.path.join(output_dir, 'compare_venn2_species_detail.csv')

        save_species_detail_table(
            {
                'only_species_list': species_list_set - mf_species_set,
                'only_MF_reference_DB': mf_species_set - species_list_set,
                'common_species': species_list_set & mf_species_set
            },
            species_detail_file
        )

        print('ploting venn2')
        draw_venn2(species_list_set, mf_species_set, output_dir)
        return mf_df, mf_unique_df

    formation = str(formation).upper()

    third_df = read_third_database(formation, format_path)

    if third_df is False:
        return

    third_unique_df = build_unique_rank_table(third_df)

    third_unique_df.to_csv(
        os.path.join(output_dir, f'{formation}_standard_taxonomy_table.csv'),
        index=False,
        encoding='utf-8-sig'
    )

    overlap_df = compare_taxonomy_overlap(mf_df, third_df)

    overlap_df.to_csv(
        os.path.join(output_dir, f'compare_{formation}_taxonomy_overlap.csv'),
        index=False,
        encoding='utf-8-sig'
    )

    third_species_set = get_rank_set(third_df, 'species')

    species_detail_file = os.path.join(output_dir, f'compare_{formation}_venn3_species_detail.csv')

    save_species_detail_table(
        {
            'only_species_list': species_list_set - mf_species_set - third_species_set,
            'only_MF_reference_DB': mf_species_set - species_list_set - third_species_set,
            f'only_{formation}_DB': third_species_set - species_list_set - mf_species_set,
            'species_list_and_MF_reference_DB': (species_list_set & mf_species_set) - third_species_set,
            f'species_list_and_{formation}_DB': (species_list_set & third_species_set) - mf_species_set,
            f'MF_reference_DB_and_{formation}_DB': (mf_species_set & third_species_set) - species_list_set,
            'common_three_datasets': species_list_set & mf_species_set & third_species_set
        },
        species_detail_file
    )

    print('ploting venn3')
    draw_venn3(species_list_set, mf_species_set, third_species_set, formation, output_dir)

    return mf_df, third_df, overlap_df, mf_unique_df, third_unique_df


