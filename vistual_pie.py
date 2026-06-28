import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import squarify
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from  matplotlib_venn import venn2
from  matplotlib_venn import venn3
from Bio import SeqIO
import re
from matplotlib.patches import Rectangle


def pie(file):
    print('ploting pie')
    if os.path.exists(file):
        try:
            csv_df = pd.read_csv(file, encoding='utf-8-sig')
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return


    # 输出
    output_fp = os.path.dirname(file)
    output_fp_csv = os.path.join(output_fp,'species_condition.csv')
    output_fp_fig = os.path.join(output_fp,'pie.png')
    # 增加图形高度，为底部文本框留出空间
    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    ax1, ax2 = axes[0], axes[1]

    # 调整子图间距，给图例和文本框留出空间
    plt.subplots_adjust(left=0.08, right=0.92, top=0.85, bottom=0.25, wspace=0.3)

    # 数据统计
    # 可扩增的序列
    df_amplify_name = []
    df_amplify_condition = []
    # 可扩增的物种
    df_amplify_species = {}
    # 物种层面
    df_concerned = set()
    set1_allcanbe = set()
    set2_cant = set()

    amplicon_length = 0

    for index, row in csv_df.iterrows():
        if row['amplify'] == 1:
            amplicon_length = int(row['amplicon_length']) + amplicon_length
            df_amplify_name.append(row['name'])
            df_amplify_condition.append(1)

            # 再判断历史情况
            if row['scientific name'] in df_amplify_species.keys():
                if df_amplify_species[row['scientific name']] == 0:
                    df_amplify_species[row['scientific name']] = 1
                    df_concerned.add(row['scientific name'])
            else:
                df_amplify_species[row['scientific name']] = 1
        else:
            # 序列不可扩增
            df_amplify_name.append(row['name'])
            df_amplify_condition.append(0)
            # 再判断物种是否可扩增
            if row['scientific name'] in df_amplify_species.keys():
                if df_amplify_species[row['scientific name']] == 1:
                    df_concerned.add(row['scientific name'])
                    continue
            else:
                df_amplify_species[row['scientific name']] = 0

    for index, row in csv_df.iterrows():
        if row['amplify'] == 1:
            if row['scientific name'] in df_concerned:
                continue
            else:
                set1_allcanbe.add(row['scientific name'])
        else:
            set2_cant.add(row['scientific name'])

    # 找出三个集合中最长的长度，以便对齐
    max_len = max(len(set1_allcanbe), len(set2_cant), len(df_concerned))

    # 把集合转成列表，不足的部分用 '' 补齐（或填 NaN）
    col1 = sorted(set1_allcanbe) + [''] * (max_len - len(set1_allcanbe))
    col2 = sorted(set2_cant) + [''] * (max_len - len(set2_cant))
    col3 = sorted(df_concerned) + [''] * (max_len - len(df_concerned))

    # 组装：表头 + 每一行
    sp_df = [['Good situation', 'Cant be amplified', 'Need concerned']]
    for i in range(max_len):
        sp_df.append([col1[i], col2[i], col3[i]])

    df = pd.DataFrame(sp_df[1:], columns=sp_df[0])
    df.to_csv(output_fp_csv)


    # 计算可扩增和不可扩增序列数量
    value1 = sum(df_amplify_condition)
    value2 = len(df_amplify_condition) - value1
    bar1 = [value1, value2]

    colors1 = ['steelblue', 'gold']
    colors2 = ['steelblue' if value == 1 else 'gold' for value in df_amplify_species.values()]

    # 绘制第一个饼图
    wedge_props1 = {'edgecolor': 'black', 'linestyle': '-', 'linewidth': 0.8, 'alpha': 0.9}
    # 绘制饼图（不显示外部标签）
    wedges1, texts1, autotexts1 = ax1.pie(
        bar1,
        labels=['', ''],
        colors=colors1,
        autopct='%1.1f%%',
        startangle=90,
        pctdistance=0.8,
        wedgeprops=wedge_props1
    )

    # 创建第一个饼图的图例 - 放在饼图内部右上角
    legend_labels1 = [f'{value1} sequences can be amplified',
                      f'{value2} sequences can not be amplified']
    ax1.legend(wedges1, legend_labels1,
               loc='upper right',
               bbox_to_anchor=(0.98, 0.98),
               frameon=True,
               fancybox=True,
               borderpad=1,
               fontsize=10)

    ax1.set_title('ISPCR Sequence Result', fontsize=16, fontweight='bold', pad=20)

    # 设置第一个饼图的百分比文本颜色
    for i, autotext in enumerate(autotexts1):
        if i == 0:  # 可扩增
            autotext.set_color('white')
        else:  # 不可扩增
            autotext.set_color('black')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(12)

    # 绘制第二个饼图
    species_values = list(df_amplify_species.values())
    species_keys = list(df_amplify_species.keys())

    value3 = int(sum(df_amplify_species.values()))
    value4 = len(species_values) - value3
    bar2 = [value3, value4]

    wedge_props2 = {'edgecolor': 'black', 'linestyle': '-', 'linewidth': 0.8, 'alpha': 0.9}
    wedges2, texts2, autotexts2 = ax2.pie(
        bar2,
        colors=colors1,
        autopct='%1.1f%%',
        startangle=90,
        wedgeprops=wedge_props2,
        pctdistance=0.7
    )
    ax2.set_title('ISPCR Species Results', fontsize=16, fontweight='bold', pad=20)

    amplify_species_num = 0
    for i in df_amplify_species.values():
        if i == 1:
            amplify_species_num += 1
    non_amplify_species_num = 0
    for i in df_amplify_species.values():
        if i == 0:
            non_amplify_species_num += 1

    # 第二个饼图的图例 - 放在饼图内部右上角
    patch1 = mpatches.Patch(color='steelblue', label=f'{amplify_species_num} species can be amplified')
    patch2 = mpatches.Patch(color='gold', label=f'{non_amplify_species_num} species can not be amplified')
    ax2.legend(handles=[patch1, patch2],
               loc='upper right',
               bbox_to_anchor=(0.98, 0.98),
               frameon=True,
               fancybox=True,
               borderpad=1,
               fontsize=10)

    # 设置第二个饼图的百分比文本颜色
    for i, (wedge, autotext) in enumerate(zip(wedges2, autotexts2)):
        if i == 0:  # 第一个扇形（可扩增物种）
            autotext.set_color('white')
        else:  # 第二个扇形（不可扩增物种）
            autotext.set_color('black')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(12)

    # 创建需要关注的物种列表
    legend_text1 = ['Species with both amplicons and non-amplicons (need attention):\n']
    df_concerned = list(df_concerned)
    if df_concerned:
        for i in range(0, min(len(df_concerned), 10)):  # 只显示前10个
            legend_text1.append(f'• {df_concerned[i]}\n')
        if len(df_concerned) > 10:
            legend_text1.append(f'... and {len(df_concerned) - 10} more species')
    else:
        legend_text1.append('None')

    legend_text1 = ' '.join(legend_text1)

    # 分类统计
    count_df = csv_df[csv_df['amplify'] == 1]
    g_num = len(count_df['genus'].value_counts().values)
    f_num = len(count_df['family'].value_counts().values)
    o_sum = len(count_df['order'].value_counts().values)
    c_sum = len(count_df['class'].value_counts().values)
    p_sum = len(count_df['phylum'].value_counts().values)
    k_sum = len(count_df['kingdom'].value_counts().values)

    # 计算平均扩增子长度
    if value1 > 0:
        avg_amplicon_length = amplicon_length / value1
    else:
        avg_amplicon_length = 0

    # 创建三个文本框 - 放在图形底部，水平排列

    # 第一个文本框：关注物种
    fig.text(0.08, 0.15, legend_text1, fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#E6F3FF", alpha=0.9, edgecolor='#1E88E5'))

    # 第二个文本框：平均扩增子长度
    length_text = f'Average Amplicon Length:\n{avg_amplicon_length:.1f} bp'
    fig.text(0.38, 0.25, length_text, fontsize=11, verticalalignment='top',
             color='red', fontweight='bold',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFF9C4", alpha=0.9, edgecolor='#FFD600'))

    # 第三个文本框：分类统计
    legend_text3 = f'Amplified Taxonomic Classification:\n' \
                   f'• Species: {amplify_species_num}\n' \
                   f'• Genus: {g_num}\n' \
                   f'• Family: {f_num}\n' \
                   f'• Order: {o_sum}\n' \
                   f'• Class: {c_sum}\n' \
                   f'• Phylum: {p_sum}\n' \
                   f'• Kingdom: {k_sum}'

    fig.text(0.68, 0.15, legend_text3, fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#E8F5E8", alpha=0.9, edgecolor='#4CAF50'))

    # 添加整体标题
    plt.suptitle('ISPCR Analysis Results', fontsize=18, fontweight='bold', y=0.95)

    # output
    #output_dir = os.path.dirname(file)
    #output_path = 'D:/1/researches/中期/2/last/dq/Cluster/pie.png'
    plt.savefig(output_fp_fig, format='png', dpi=800, bbox_inches='tight')
    print(f"pie figure saved in: {output_fp_fig}")
    print(f'species condition saved in: {output_fp_csv}')
    print('\n')
    print('Sequences to filter or not NEED combine Cluster results')


pie(file=r'D:\1\researches\中期\DJ\dj\Filter_T\Filter.csv')