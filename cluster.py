
from Bio import SeqIO
import os
import pandas as pd
import subprocess



def cluster_sequences(input_fasta, threshold):
    uc_file = "clusters.uc"
    vsearch_path = r'F:\gdut\学习\2中期\2\vsearch-2.30.0-win-x86_64\bin\vsearch.exe'
    cmd = [
        vsearch_path,
        "--cluster_size", input_fasta,
        "--id", str(threshold),
        "--uc", uc_file,
        "--threads", str(32),
        "--maxseqlength", str(999999)
    ]
    try:
        subprocess.run(cmd, check=True)
        # print("Vsearch is finished")
    except subprocess.CalledProcessError as e:
        print(f"cluster fall: {e}")
        return []

    #数据过大
    if not os.path.exists(uc_file) or os.path.getsize(uc_file) == 0:
        cluster_info = []
        for seq_record in SeqIO.parse(input_fasta,format='fasta'):
            cluster_info.append({
                "name": seq_record.id,
                "represent": f"na"
            })
        return cluster_info



    file = pd.read_table(uc_file, encoding='utf-8-sig', header=None, sep='\t')
    # print(file)
    cluster_info = []
    cluster_name = {}
    for index, row in file.iterrows():
        if row[0] == 'C':
            cluster_name[row[8]] = int(row[2])
    sorted_by_value_asc = sorted(cluster_name.items(), key=lambda x: x[1], reverse=True)
    # print(sorted_by_value_asc)
    sorted_dict = dict(sorted_by_value_asc)

    i = 1
    for key in sorted_dict.keys():

        cluster_info.append({
            "name": key,
            "represent": f"center_{i}"
        })
        for index, row in file.iterrows():
            if row[0] == 'H' and row[9] == key:
                cluster_info.append({
                    "name": row[8],
                    "represent": f"clustered_{i}"
                })
        i = i + 1

    return cluster_info


def cluster(csv_path, file_path, threshold):
    # system path
    if os.path.exists(csv_path) and os.path.exists(file_path) == True:
        df = pd.read_csv(csv_path)
        # get same path
        print('please put csv and file in the same second level directory, and cluster file will output in the same')
        if os.path.dirname(os.path.dirname(csv_path)) == os.path.dirname(os.path.dirname(file_path)):
            output_dir = os.path.dirname(os.path.dirname(csv_path))
            print('read csv and file successfully')

        else:
            print('csv and file are not in the same ')
            return

        output_dir_file = os.path.join(output_dir, 'Cluster')
        os.makedirs(output_dir_file, exist_ok=True)

        cluster_info = []
        combined_seqs = []

        files = os.listdir(file_path)
        # print(files)
        for file in files:
            input_path_cluster = f'{file_path}/{file}'
            if file.endswith('.fasta') == False:
                print(f'{file_path} has wrong {file} file，ignored')
            else:
                cluster_info.extend(cluster_sequences(input_fasta=input_path_cluster,threshold=threshold))
                for seq_record in SeqIO.parse(input_path_cluster, format='fasta'):
                    combined_seqs.append(seq_record)



        df_cluster_info = pd.DataFrame(cluster_info)
        #print(df_cluster_info)


        df_merge = pd.merge(df,df_cluster_info,on='name',how='left')

        for index, row in df_merge.iterrows():
            if row['represent'] == '':
                df_merge.loc[index,'represent'] = uncluster_999

        output_dir_cluster_csv = os.path.join(output_dir_file,'Cluster_info.csv')
        df_merge.to_csv(output_dir_cluster_csv, index=False , encoding="utf-8-sig")

        output_dir_cluster_seq = os.path.join(output_dir_file, 'Combined_Cluster.fasta')
        SeqIO.write(combined_seqs, output_dir_cluster_seq, format='fasta')


        #添加无法过滤的检查







    else:
        print('files do not exist')
        return


cluster(file_path=r'F:\gdut\学习\3答辩\嵌入\dq\ISPCR\ISPCR_sequence',
        csv_path=r'F:\gdut\学习\3答辩\嵌入\dq\ISPCR\ISPCR_info.csv', threshold=0.97)


