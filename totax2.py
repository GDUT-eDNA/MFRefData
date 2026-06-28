import pandas as pd
import requests
import os
import tarfile
#
def download_NCBI_dump(path):

    return print('download dump files successfully')

def totax(csv, dump, download_dump = False,download_path = None):
    if download_dump and download_path == True:
        download_NCBI_dump(download_path)
        return

    if csv:
        if os.path.exists(csv) == False:
            print('pleasure ensure csv files exist')
            return

        if  os.path.exists(dump) == True:
            names_path = os.path.join(dump, 'names.dmp')
            nodes_path = os.path.join(dump, 'nodes.dmp')
            if os.path.exists(names_path) & os.path.exists(nodes_path) == True:
                print('prepare reload dump')
            else:
                return print('dump files is missing')
            input_csv = pd.read_csv(csv)
            #
            output_path_dir = os.path.dirname(os.path.dirname(csv))
            output_path = os.path.join(output_path_dir, 'Totax')
            os.makedirs(output_path, exist_ok= True)
            #
            print('this step is compatible with NCBI sequences with taxid')
            print('sequences without taxid value will be junmped')
            print('\n')
            print('reload csv file successfully')
            name_dmp = pd.read_table(names_path, sep='\t\\|\t',engine='python', header=None)
            print('reload names.dmp file successful')
            node_dmp = pd.read_table(nodes_path, sep='\t\\|\t',engine='python', header=None)
            print('reload nodes.dmp file successful')

            # 修改点1: 确保node_dmp的第一列是数值类型，以便正确匹配
            node_dmp[0] = pd.to_numeric(node_dmp[0], errors='coerce')
            node_dmp[1] = pd.to_numeric(node_dmp[1], errors='coerce')

            #建立集合
            list_lineage_all = []

            # 去除重复taxid
            taxid_dereplicate = set()
            for taxid in input_csv['taxid']:
                taxid_dereplicate.add(int(taxid))  # 修改点2: 确保taxid为整数

            # 去重复后的taxid查找父类节点full lineage
            for taxid in taxid_dereplicate:
                list_lineage_single = []
                if taxid not in node_dmp[0].values:  # 修改点3: 使用.values进行正确检查
                    print(f'{taxid}不在该nodes文件里，请更新tax文件或采取其他方式查询注释信息')
                else:
                    current_taxid = taxid  # 修改点4: 使用临时变量，避免修改循环变量
                    while current_taxid != 1:
                        #找到行
                        tax_row = node_dmp[node_dmp[0] == current_taxid]
                        # 修改点5: 添加空值检查
                        if tax_row.empty:
                            print(f'Warning: taxid {current_taxid} not found during lineage construction')
                            break
                        tax_level = tax_row.iloc[0, 2]
                        list_lineage_single.append((current_taxid, tax_level))
                        current_taxid = int(tax_row.iloc[0, 1])
                    # 修改点6: 只有在成功到达根节点时才添加root
                    if current_taxid == 1:
                        list_lineage_single.append((1, 'all_root'))
                list_lineage_all.append(list_lineage_single[::-1])

            # for i in list_lineage_all:
            #     print(i)

            # 在所有的节点中去除重复nodes
            taxid_sets = set()
            # 构建字典，方便查询taxid
            tax_to_name = {}

            for i in list_lineage_all:
                for j in i:
                    taxid_in_single_node = j[0]
                    taxid_sets.add(taxid_in_single_node)

            # 用除去重复的nodes查询name,以字典储存每个节点的taxid和name
            for taxid_set in taxid_sets:
                name_row = name_dmp[(name_dmp[0] == taxid_set)]
                if (name_row[3] == str('scientific name\t|')).any():
                    name_row = name_row[name_row[3] == str('scientific name\t|')]
                    # print(name_row)#注意，提取的是series格式还是iloc的值格式
                    # 修改点7: 添加空值检查
                    if not name_row.empty:
                        tax_to_name[taxid_set] = name_row.iloc[0, 1]
                    else:
                        tax_to_name[taxid_set] = 'Unknown'

            # print(tax_to_name)
            # 使用新列表存储，元组不可改变
            list_lineage_all_with_name = []

            # 匹配结果，添加名称
            for lineage in list_lineage_all:
                list_lineage_with_name = []
                for j in lineage:
                    taxid = j[0]
                    tax_level = j[1]
                    if taxid in tax_to_name:
                        new_lineage = (taxid, tax_level, tax_to_name[taxid])
                    else:
                        new_lineage = (taxid, tax_level, 'NA')
                    list_lineage_with_name.append(new_lineage)
                # print(list_lineage_with_name)
                list_lineage_all_with_name.append(list_lineage_with_name)


            # for i in list_lineage_all_with_name:
            #     print(i)

            #输出匹配full_lineage
            dict_taxid_lineage = {}
            for i, lineage in enumerate(list_lineage_all_with_name):
                # 获取对应的taxid（按顺序）
                taxid = list(taxid_dereplicate)[i]
                # 格式化lineage
                formatted_lineage = []
                for node in lineage:
                    formatted_lineage.append(f'{node[0]}_{node[1]}_{node[2]}')
                #使用空格连接所有节点
                dict_taxid_lineage[taxid] = lineage#'  '.join(formatted_lineage)


            # 输出指定lineage
            input_csv['kingdom'] = ''
            input_csv['phylum'] = ''
            input_csv['class'] = ''
            input_csv['order'] = ''
            input_csv['family'] = ''
            input_csv['genus'] = ''
            input_csv['species'] = ''

            input_csv['full_lineage'] = ''


            for index, row in input_csv.iterrows():
                taxid = row['taxid']
                #print(taxid)
                for lineage in list_lineage_all_with_name:
                    # 先判断是否为空，因为lineage有空值
                    if not lineage:
                        print(f'{lineage} is empty')
                        continue

                    # 根据taxid查找结果的最后一项都为taxid，再在有数据的里面进行匹配
                    if taxid != lineage[-1][0]:
                        continue
                        #print(f'species {taxid} is missing in {lineage}')

                    if taxid == lineage[-1][0]:
                        input_csv.at[index, 'kingdom'] = 'na'
                        input_csv.at[index, 'phylum'] = 'na'
                        input_csv.at[index, 'class'] = 'na'
                        input_csv.at[index, 'order'] = 'na'
                        input_csv.at[index, 'family'] = 'na'
                        input_csv.at[index, 'genus'] = 'na'
                        input_csv.at[index, 'species'] = 'na'

                        full_lineage_fomatted = []
                        for node in lineage:
                            if node[1] == 'kingdom':
                                input_csv.at[index, 'kingdom'] = (f'k_{node[2]}')
                            if node[1] == 'phylum':
                                input_csv.at[index, 'phylum'] = (f'p_{node[2]}')
                            if node[1] == 'class':
                                input_csv.at[index, 'class'] = (f'c_{node[2]}')
                            if node[1] == 'order':
                                input_csv.at[index, 'order'] = (f'o_{node[2]}')
                            if node[1] == 'family':
                                input_csv.at[index, 'family'] = (f'f_{node[2]}')
                            if node[1] == 'genus':
                                input_csv.at[index, 'genus'] = (f'g_{node[2]}')
                            if node[1] == 'species':
                                input_csv.at[index, 'species'] = (f's_{node[2]}')
                            # full_linage
                            full_lineage_fomatted.append((f'{node[0]}_{node[1]}_{node[2]};'))
                        full_lineage_fomatted = ''.join(full_lineage_fomatted)
                        # print(full_lineage_fomatted)
                        # print(type(full_lineage_fomatted))
                        input_csv.at[index, 'full_lineage'] = full_lineage_fomatted


            for index, row in input_csv.iterrows():
                if not row['class'] and not row['kingdom'] and not row['phylum'] and not row['family'] and not row['species'] and not row['genus'] and not row['order']:
                    input_csv.loc[index,'kingdom'] = 'na'
                    input_csv.loc[index,'phylum'] = 'na'
                    input_csv.loc[index,'class'] = 'na'
                    input_csv.loc[index,'order'] = 'na'
                    input_csv.loc[index,'family'] = 'na'
                    input_csv.loc[index,'genus'] = 'na'
                    input_csv.loc[index,'species'] = 'na'
                    input_csv.loc[index, 'full_lineage'] = 'na'



    input_csv.to_csv(f'{output_path}/Totax.csv', index=False, encoding='utf_8_sig')
    print('finish totax')

totax(csv = r'F:\gdut\学习\3答辩\嵌入\dq\Cluster\Cluster_info.csv', dump= r'F:\gdut\学习\2中期\1\taxdmp')