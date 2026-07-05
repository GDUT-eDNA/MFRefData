import pandas as pd
from Bio import SeqIO
import os
from pathlib import Path


# 可实现多参数过滤
# 过滤逻辑
# 基于序列质量得分过滤
# 基于指定行列过滤
# 基于序列长度过滤
# 基于聚类情况
# 1.过滤所有序列中心，但保留序列簇
# 2.过滤所有序列簇，但保留中心
# 3.过滤所有指定level以下的序列中心和簇
# 基于某个序列注释信息进行过滤
def filter(csv_file, fasta_file,
           column=False, value=False,
           tax_miss=False, low_qual=False,
           no_center=False, no_cluster=False,
           no_cluster_level=False,
           min_len=False, max_len=False,
           raw_max_len=False,
           keep_level=False):

    # --------------------------------------------------
    # 修改部分：分别检查 CSV 和 FASTA 文件
    # --------------------------------------------------
    if not os.path.isfile(csv_file):
        print(f'csv filepath is wrong: {csv_file}')
        return pd.DataFrame

    if not os.path.isfile(fasta_file):
        print(f'fasta filepath is wrong: {fasta_file}')
        return pd.DataFrame

    if not csv_file.endswith('.csv'):
        print('csv suffix is wrong, please check')
        return pd.DataFrame

    if not fasta_file.endswith('.fasta'):
        print('fasta suffix is wrong, please check')
        return pd.DataFrame

    # 保留原变量名称，避免修改后续 FASTA 逻辑
    sequences_file = fasta_file

    try:
        csv_data = pd.read_csv(csv_file)
        print('read csv file successfully ')
    except Exception as e:
        print('csv filepath is wrong')
        return pd.DataFrame

    # 总体过滤列表
    sequences_toremove = set()
    out_put_path_dir = {}

    # --------------------------------------------------
    # 低质量过滤
    # --------------------------------------------------
    if low_qual != False:
        print('caculate low quality')
        out_put_path_dir['low_qual'] = low_qual

        df = csv_data[
            [
                'taxid',
                'name',
                'represent',
                'amplify',
                'geography',
                'reference'
            ]
        ]

        taxid_list = csv_data[
            'taxid'
        ].value_counts().index.tolist()

        for taxid in taxid_list:
            small_df = df[df['taxid'] == taxid]

            for index, row in small_df.iterrows():
                score = 0

                if row['geography'] != 'na':
                    score = score + 1

                represent_info = row['represent']
                represent_info_count = 0

                for index2, row2 in small_df.iterrows():
                    if row2['represent'] == represent_info:
                        represent_info_count = (
                            represent_info_count + 1
                        )

                score = (
                    score
                    + int(row['reference'])
                    + represent_info_count
                    + int(row['amplify'] * 2)
                )

                if score < low_qual:
                    sequences_toremove.add(
                        str(row['name'])
                    )
                else:
                    continue

    else:
        out_put_path_dir['low_qual'] = '_'

    # --------------------------------------------------
    # 分类信息缺失过滤
    # --------------------------------------------------
    if tax_miss != False:
        print('caculate NA')
        out_put_path_dir['tax_miss'] = tax_miss

        df = csv_data[
            [
                'taxid',
                'name',
                'kingdom',
                'phylum',
                'class',
                'order',
                'family',
                'genus',
                'species'
            ]
        ]

        list_concern = [
            'kingdom',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species'
        ]

        for index, row in df.iterrows():
            miss_count = 0

            for head in list_concern:
                value_ = row[head]

                if value_ == 'na':
                    miss_count += 1
                else:
                    continue

            if miss_count > tax_miss:
                sequences_toremove.add(row['name'])

    else:
        out_put_path_dir['tax_miss'] = '_'

    # --------------------------------------------------
    # 指定列、指定值过滤
    # --------------------------------------------------
    if column and value:
        out_put_path_dir['column_value'] = f'{column}_'

        if isinstance(value, str) and value.endswith('.txt'):
            if os.path.exists(value) == True:
                toremove = pd.read_csv(
                    value,
                    sep="\t",
                    header=None,
                    encoding='utf-8'
                )[0].tolist()

                for i in toremove:
                    sequences_toremove.update(
                        map(
                            str,
                            csv_data[
                                csv_data[column] == str(i)
                            ]['name'].values
                        )
                    )

        else:
            sequences_toremove.update(
                map(
                    str,
                    csv_data[
                        csv_data[column] == value
                    ]['name'].values
                )
            )

    else:
        out_put_path_dir['column_value'] = '_'

    # --------------------------------------------------
    # 保留指定分类等级
    # --------------------------------------------------
    if keep_level:
        out_put_path_dir['keep_level'] = keep_level

        df = csv_data[
            [
                'taxid',
                'name',
                'kingdom',
                'phylum',
                'class',
                'order',
                'family',
                'genus',
                'species'
            ]
        ]

        for index, row in df.iterrows():
            if row[keep_level] == 'na':
                sequences_toremove.add(row['name'])

    else:
        out_put_path_dir['keep_level'] = '_'

    # --------------------------------------------------
    # 只保留簇序列：删除中心序列
    # --------------------------------------------------
    if no_center:
        out_put_path_dir['no_center'] = 'T'

        for index, row in csv_data.iterrows():
            if (
                str(row['represent']).split('_')[0]
                == 'center'
            ):
                sequences_toremove.add(row['name'])
            else:
                continue

    else:
        out_put_path_dir['no_center'] = '_'

    # --------------------------------------------------
    # 只保留中心：删除簇序列
    # --------------------------------------------------
    if no_cluster:
        # 保留原代码写法
        out_put_path_dir['no_center'] = 'T'

        for index, row in csv_data.iterrows():
            if (
                str(row['represent']).split('_')[0]
                == 'cluster'
            ):
                sequences_toremove.add(row['name'])
            else:
                continue

    else:
        out_put_path_dir['no_cluster'] = '_'

    # --------------------------------------------------
    # 聚类等级过滤
    # --------------------------------------------------
    if no_cluster_level:
        # 保留原代码写法
        out_put_path_dir['no_center'] = 'T'

        for index, row in csv_data.iterrows():
            if (
                str(row['represent']).split('_')[1]
                > str(no_cluster_level)
            ):
                sequences_toremove.add(row['name'])
            else:
                continue

    else:
        out_put_path_dir['no_cluster_level'] = '_'

    # --------------------------------------------------
    # 最小扩增子长度过滤
    # --------------------------------------------------
    if min_len:
        # 保留原代码写法
        out_put_path_dir['no_center'] = f'{min_len}'

        for index, row in csv_data.iterrows():
            temp = int(row['amplicon_length'])

            if temp != None:
                if row['amplicon_length'] < min_len:
                    sequences_toremove.add(
                        row['name']
                    )
                else:
                    continue
            else:
                return print('missing amplicon_length')

    else:
        out_put_path_dir['min_len'] = '_'

    # --------------------------------------------------
    # 最大扩增子长度过滤
    # --------------------------------------------------
    if max_len:
        # 保留原代码写法
        out_put_path_dir['no_center'] = f'{max_len}'

        for index, row in csv_data.iterrows():
            temp = int(row['amplicon_length'])

            if temp != None:
                if row['amplicon_length'] > max_len:
                    sequences_toremove.add(
                        row['name']
                    )
                else:
                    continue
            else:
                return print('missing amplicon_length')

    else:
        out_put_path_dir['max_len'] = '_'

    # --------------------------------------------------
    # 原始序列最大长度过滤
    # --------------------------------------------------
    if raw_max_len:
        # 保留原代码写法
        out_put_path_dir['no_center'] = (
            f'{raw_max_len}'
        )

        for index, row in csv_data.iterrows():
            if row['raw_seq_length']:
                if row['raw_seq_length'] > raw_max_len:
                    sequences_toremove.add(
                        row['name']
                    )
                else:
                    continue
            else:
                return print('missing raw_seq_length')

    else:
        out_put_path_dir['raw_max_len'] = '_'

    # --------------------------------------------------
    # 生成名称
    # --------------------------------------------------
    output_name = 'Filter'

    for key, value in out_put_path_dir.items():
        if value != '_':
            output_name = (
                output_name + f'_{value}'
            )
        else:
            continue

    print(output_name)

    # --------------------------------------------------
    # 输出路径
    # --------------------------------------------------
    input_file_dir = os.path.dirname(
        os.path.abspath(csv_file)
    )

    output_path = os.path.dirname(
        input_file_dir
    )

    out_put_path = os.path.join(
        output_path,
        str(output_name)
    )

    output_dir_seq = os.path.join(
        out_put_path,
        'Filter.fasta'
    )

    output_dir_csv = os.path.join(
        out_put_path,
        'Filter.csv'
    )

    os.makedirs(
        out_put_path,
        exist_ok=True
    )

    # --------------------------------------------------
    # 读取和过滤 FASTA
    # --------------------------------------------------
    new_seq_record = []

    for seq_record in SeqIO.parse(
        sequences_file,
        format='fasta'
    ):
        if seq_record.id in sequences_toremove:
            continue
        else:
            new_seq_record.append(seq_record)

    SeqIO.write(
        new_seq_record,
        output_dir_seq,
        'fasta'
    )

    # --------------------------------------------------
    # 过滤和输出 CSV
    # --------------------------------------------------
    new_df = csv_data[
        ~csv_data['name'].isin(
            sequences_toremove
        )
    ]

    print(
        f'filtering '
        f'{len(sequences_toremove)} sequences'
    )

    new_df.to_csv(
        output_dir_csv,
        encoding='utf-8',
        index=False
    )


if __name__ == "__main__":
    work_dir = './work'
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    filter(
        csv_file=os.path.join(work_dir, 'Totax', 'Totax.csv'),
        fasta_file=os.path.join(work_dir, 'Cluster', 'Combined_Cluster.fasta'),
        raw_max_len=5000
    )
