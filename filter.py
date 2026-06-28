
import pandas as pd
from Bio import SeqIO
import os

###可实现多参数过滤
#过滤逻辑
#基于序列质量得分过滤
#基于指定行列过滤
#基于序列长度过滤
#基于聚类情况
    #1.过滤所有序列中心，但保留序列簇
    #2.过滤所有序列簇，但保留中心
    #3.过滤所有指定level以下的序列中心和簇
#基于某个序列注释信息进行过滤
def filter(file, column = False, value = False, tax_miss = False, low_qual = False,
           no_center = False,no_cluster = False, no_cluster_level = False,
           min_len = False, max_len = False,raw_max_len = False,
           keep_level = False):

    #打开文件
    if os.path.exists(file) == True:
        csv_file = None
        sequences_file = None
        for filename in os.listdir(file):
            file_path = os.path.join(file,filename)
            if file_path.endswith('.csv'):
                csv_file = file_path

            elif file_path.endswith('.fasta'):
                sequences_file = file_path
        if csv_file is None or sequences_file is None:
            print('csv or fasta is wrong, please check')


        try:
            csv_data = pd.read_csv(csv_file)
            print('read csv file successfully ')
        except Exception as e:
            print('csv filepath is wrong')
            return pd.DataFrame

        #总体过滤列表
        sequences_toremove = set()
        out_put_path_dir = {}



        ###低质量过滤
        if low_qual != False:
            print('caculate low quality')
            out_put_path_dir['low_qual'] = low_qual


            df = csv_data[['taxid', 'name', 'represent', 'amplify', 'geography', 'reference']]
            taxid_list = csv_data['taxid'].value_counts().index.tolist()


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
                            represent_info_count = represent_info_count + 1
                    score = score + int(row['reference']) + represent_info_count + int(row['amplify']*2)
                    if score < low_qual:
                        sequences_toremove.add(str(row['name']))
                    else:
                        continue

        else:
            out_put_path_dir['low_qual'] = '_'


        ###
        if tax_miss != False:
            print('caculate NA')
            out_put_path_dir['tax_miss'] = tax_miss


            df = csv_data[['taxid', 'name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]
            list_concern = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            for index,row in df.iterrows():
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



        ###
        if column and value:
            out_put_path_dir['column_value'] = f'{column}_'
            if isinstance(value, str) and value.endswith('.txt'):
                if os.path.exists(value) == True:
                    toremove = pd.read_csv(value, sep="\t", header=None, encoding='utf-8')[0].tolist()
                    for i in toremove:
                        sequences_toremove.update(map(str, csv_data[csv_data[column] == str(i)]['name'].values))

            else:
                sequences_toremove.update(map(str, csv_data[csv_data[column] == value]['name'].values))

        else:
            out_put_path_dir['column_value'] = '_'


        ###
        if keep_level:
            out_put_path_dir['keep_level'] = keep_level

            df = csv_data[['taxid', 'name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]
            for index, row in df.iterrows():
                if row[keep_level] == 'na':
                    sequences_toremove.add(row['name'])
        else:
            out_put_path_dir['keep_level'] = '_'


        ###只保留簇序列
        if no_center:
            out_put_path_dir['no_center'] = 'T'
            for index, row in csv_data.iterrows():
                if str(row['represent']).split('_')[0] == 'center':
                    sequences_toremove.add(row['name'])
                else:
                    continue
        else:
            out_put_path_dir['no_center'] = '_'

        ###只保留中心
        if no_cluster:
            out_put_path_dir['no_center'] = 'T'
            for index, row in csv_data.iterrows():
                if str(row['represent']).split('_')[0] == 'cluster':
                    sequences_toremove.add(row['name'])
                else:
                    continue
        else:
            out_put_path_dir['no_cluster'] = '_'

        ###保留等级
        if no_cluster_level:
            out_put_path_dir['no_center'] = 'T'
            for index, row in csv_data.iterrows():
                if str(row['represent']).split('_')[1] > str(no_cluster_level):
                    sequences_toremove.add(row['name'])
                else:
                    continue
        else:
            out_put_path_dir['no_cluster_level'] = '_'

        ###保留一定长度
        if min_len:
            out_put_path_dir['no_center'] = f'{min_len}'
            for index, row in csv_data.iterrows():
                temp = int(row['amplicon_length'])
                if temp != None:
                    if row['amplicon_length'] < min_len:
                        sequences_toremove.add(row['name'])
                    else:
                        continue
                else:
                    return print('missing amplicon_length')
        else:
            out_put_path_dir['min_len'] = '_'

        if max_len:
            out_put_path_dir['no_center'] = f'{max_len}'
            for index, row in csv_data.iterrows():
                temp = int(row['amplicon_length'])
                if temp != None:
                    if row['amplicon_length'] > max_len:
                        sequences_toremove.add(row['name'])
                    else:
                        continue
                else:
                    return print('missing amplicon_length')
        else:
            out_put_path_dir['max_len'] = '_'

        if raw_max_len:
            out_put_path_dir['no_center'] = f'{raw_max_len}'
            for index, row in csv_data.iterrows():
                if row['raw_seq_length']:
                    if row['raw_seq_length'] > raw_max_len:
                        sequences_toremove.add(row['name'])
                    else:
                        continue
                else:
                    return print('missing raw_seq_length')
        else:
            out_put_path_dir['raw_max_len'] = '_'


        #生成名称
        output_name = 'Filter'
        for key, value in out_put_path_dir.items():
            if value != '_':
                output_name = output_name + f'_{value}'
            else:
                continue
        print(output_name)

        output_path = os.path.dirname(file)
        out_put_path = os.path.join(output_path, str(output_name))
        output_dir_seq = os.path.join(out_put_path,'Filter.fasta')
        output_dir_csv = os.path.join(out_put_path,'Filter.csv')
        os.makedirs(out_put_path,exist_ok=True)


        new_seq_record = []

        for seq_record in SeqIO.parse(sequences_file, format='fasta'):
            if seq_record.id in sequences_toremove:
                continue
            else:
                new_seq_record.append(seq_record)
        SeqIO.write(new_seq_record, output_dir_seq, 'fasta')

        new_df = csv_data[~csv_data['name'].isin(sequences_toremove)]
        print(f'filtering {len(sequences_toremove)} sequences')
        new_df.to_csv(output_dir_csv, encoding='utf-8',index=False)




#filter(file=r'C:\Users\25327\Desktop\mock community\sequences\mf\TELE02\Cluster',raw_max_len=90000,tax_miss=2,max_len=500)
filter(file=r'C:\Users\25327\Desktop\mock community\sequences\mf\TELE02\Filter_2_90000',no_cluster_level=5)

#filter(file='D:/1/researches/中期/DJ/dj/Filter_2',no_cluster_level=5)

#filter(file='D:/1/researches/model/embed_clean/Cluster',column='name',value='C:/Users/25327/Desktop/problematic_headers.txt')
#filter(file='D:/1/researches/中期/DJ/dj_extend/Filter_scientific name_',column='genus',value='D:/1/researches/中期/DJ/remove.txt')

#filter(file=r'D:\1\researches\中期\test\Cluster',raw_max_len=5000)