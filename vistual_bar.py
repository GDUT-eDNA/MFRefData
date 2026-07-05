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
from matplotlib_venn.layout.venn2 import DefaultLayoutAlgorithm

def get_colors(n):
    colors = []
    for i in range(n):
        if i < 20:
            colors.append(plt.cm.tab20(i % 20))
        elif i < 40:
            colors.append(plt.cm.tab20c((i - 20) % 20))
        elif i < 60:
            colors.append(plt.cm.Set3((i - 40) % 12))
        else:
            # 使用HSV，黄金比例分布色相
            hue = (i * 0.618033988749895) % 1.0
            colors.append(plt.cm.hsv(hue))
    return colors


#堆叠图
def bar(inutfile,output_file):
    print('ploting')
    if os.path.exists(inutfile):
        try:
            csv_df = pd.read_csv(inutfile, encoding='utf-8-sig')
        except Exception as e:
            return


    # 创建2x2的子图画布
    fig, axes = plt.subplots(2, 2, figsize=(23, 15),constrained_layout=True)
    ax1, ax2, ax4, ax3 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # 设置全局字体和样式
    # plt.rcParams['font.size'] = 8

    # 图表1：分类组成堆叠条形图
    taxonomic_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomic_levels.reverse()


    # 获取数据
    appendix_row = [0] * 7
    remained_num = [0] * 7


    species_row = []
    species_row_count = csv_df['species'].value_counts()
    if len(species_row_count) >= 21:
        species_row.extend(species_row_count.head(21).tolist())
        appendix_row[0] = species_row_count.values[21:].sum()
        remained_num[0] = len(species_row_count) - 21
    else:
        species_row.extend(species_row_count.values.tolist())
    matrix_len = len(species_row)


    genus_row = []
    genus_row_count = csv_df['genus'].value_counts()
    if len(genus_row_count) >= matrix_len:
        genus_row.extend(genus_row_count.head(matrix_len).tolist())
        appendix_row[1] = genus_row_count.values[matrix_len:].sum()
        remained_num[1] = len(genus_row_count) - matrix_len
    elif len(genus_row_count) < matrix_len:
        genus_row.extend(genus_row_count.values.tolist())
        while len(genus_row) < matrix_len:
            genus_row.append(0)


    family_row = []
    family_row_count = csv_df['family'].value_counts()
    if len(family_row_count) >= matrix_len:
        family_row.extend(family_row_count.head(matrix_len).tolist())
        appendix_row[2] = family_row_count.values[matrix_len:].sum()
        remained_num[2] = len(family_row_count) - matrix_len
    elif len(family_row_count) < matrix_len:
        family_row.extend(family_row_count.values.tolist())
        while len(family_row) < matrix_len:
            family_row.append(0)

    order_row = []
    order_row_count = csv_df['order'].value_counts()
    if len(order_row_count) >= matrix_len:
        order_row.extend(order_row_count.head(matrix_len).tolist())
        appendix_row[3] = order_row_count.values[matrix_len:].sum()
        remained_num[3] = len(order_row_count) - matrix_len
    elif len(order_row_count) < matrix_len:
        order_row.extend(order_row_count.values.tolist())
        while len(order_row) < matrix_len:
            order_row.append(0)

    class_row = []
    class_row_count = csv_df['class'].value_counts()
    if len(class_row_count) >= matrix_len:
        class_row.extend(class_row_count.head(matrix_len).tolist())
        appendix_row[4] = class_row_count.values[matrix_len:].sum()
        remained_num[4] = len(class_row_count) - matrix_len
    elif len(class_row_count) < matrix_len:
        class_row.extend(class_row_count.values.tolist())
        while len(class_row) < matrix_len:
            class_row.append(0)

    phylum_row = []
    phylum_row_count = csv_df['phylum'].value_counts()
    if len(phylum_row_count) >= matrix_len:
        phylum_row.extend(phylum_row_count.head(matrix_len).tolist())
        appendix_row[5] = phylum_row_count.values[matrix_len:].sum()
        remained_num[5] = len(phylum_row_count) - matrix_len
    elif len(phylum_row_count) < matrix_len:
        phylum_row.extend(phylum_row_count.values.tolist())
        while len(phylum_row) < matrix_len:
            phylum_row.append(0)

    kingdom_row = []
    kingdom_row_count = csv_df['kingdom'].value_counts()
    if len(kingdom_row_count) >= matrix_len:
        kingdom_row.extend(kingdom_row_count.head(matrix_len).tolist())
        appendix_row[6] = kingdom_row_count.values[matrix_len:].sum()
        remained_num[6] = len(kingdom_row_count) - matrix_len
    elif len(kingdom_row_count) < matrix_len:
        kingdom_row.extend(kingdom_row_count.values.tolist())
        while len(kingdom_row) < matrix_len:
            kingdom_row.append(0)


    matrix_df = [species_row, genus_row, family_row, order_row, class_row, phylum_row, kingdom_row]
    #print(matrix_df)

    #图1的颜色
    colors = plt.cm.tab20(np.linspace(0, 1, matrix_len))
    bottom = np.zeros(len(taxonomic_levels))
    matrix_array = np.array(matrix_df)

    # 设置图表1边框 - 只保留XY轴
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)

    ax1.set_title(f'Sequences Composition(Colored Top {matrix_len})',fontsize = 15)
    ax1.set_xlabel('Taxonomic Level',fontsize = 12)
    ax1.set_ylabel('Number of Sequence',fontsize = 12)

    # 存储每个堆叠部分的高度和位置，用于添加数量标签
    bar_heights = np.zeros((len(taxonomic_levels), matrix_len))

    #一行一行堆上去，i为堆的次数
    #记录最后一行
    bottom_last = [[]]
    for i in range(matrix_len):
        #前面是：行整体切片，i为列（次数）
        row = matrix_array[:, i]
        bars = ax1.bar(taxonomic_levels, row, bottom=bottom,
                       color=colors[i], alpha=0.9, edgecolor='black', linewidth=0.3)
        # 记录每个堆叠部分的高度和位置
        for j, bar in enumerate(bars):
            bar_heights[j, i] = row[j]
        # 堆叠的原理
        bottom += row
        bottom_last = bottom


    # 为每个小柱子添加数量标签
    for i in range(len(taxonomic_levels)):
        total_height = bottom_last[i]
        cumulative_height = 0
        for j in range(matrix_len):
            height = bar_heights[i, j]
            if height > 0 and height / total_height > 0.009:  # 只添加非零值的标签
                # 计算标签位置（在柱子的中间）
                label_y = cumulative_height + height / 2
                ax1.text(i, label_y, f'{int(height)}',
                         ha='center', va='center', fontsize=6, fontweight='bold',
                         color='black')
            cumulative_height += height

    # 继续第二次堆叠
    bars = ax1.bar(
        taxonomic_levels, appendix_row, bottom=bottom_last,
        color = 'whitesmoke',alpha=0.8,edgecolor='black', linewidth=0.3
    )
    # bar_heights_last = bars.get_height()

    # 添加侧面标签
    for i in range(len(taxonomic_levels)):
        if appendix_row[i] == 0:
            continue
        else:
            label_y = bottom_last[i] + appendix_row[i] / 2
            ax1.text(i, label_y, f'remaining {int(appendix_row[i])} seqs\nbelonged to {remained_num[i]} {taxonomic_levels[i]}',
                         ha='center', va='center', fontsize=6, fontweight='bold',
                         color='maroon')


### 图表2：目-属堆叠条形图,等测试结果出来，看是否也需要截取
    order_name_label = []
    order_df = csv_df[['name', 'order', 'genus', 'species']]
    order_count = order_df['order'].value_counts()
    order_name_label.extend(order_count.index.tolist())
    genus_name = []
    genus_name_value = csv_df['genus'].value_counts()
    genus_name.extend(genus_name_value.index.tolist())

    order_matrix = []
    for single_genus in genus_name:
        single_row = []
        for single_order in order_name_label:
            if (csv_df[csv_df['genus'] == single_genus]['order'] == single_order).any() == True:
                single_row.append(csv_df[csv_df['genus'] == single_genus]['species'].nunique())
            else:
                single_row.append(0)
        order_matrix.append(single_row)

    # print(order_name_label)
    # print(order_matrix)
    # print(genus_name)

    matrix_array = np.array(order_matrix).T
    #print(matrix_array)




    # 使用确保区分度的函数
    colors = get_colors(len(genus_name))
    bottom = np.zeros(len(order_name_label))

    # 设置图表2边框 - 只保留XY轴
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(2)
    ax2.spines['bottom'].set_linewidth(2)


    for i, genus in enumerate(genus_name):
        counts = matrix_array[:, i]
        ax2.bar(order_name_label, counts, bottom=bottom,
                label=genus, color=colors[i], alpha=0.9, edgecolor='black', linewidth=0.5)
        bottom += counts

    ax2.set_title('Species Distribution by Order and Genus',fontsize = 15)
    ax2.set_xlabel('Order',fontsize = 12)
    ax2.set_ylabel('Number of Genus',fontsize = 12)
    ax2.set_xticklabels(order_name_label,rotation=90, ha='right', fontsize=7)

##截取标签
    legend_labels = genus_name[:len(genus_name) // 7]  # 或其他截取比例

    # 创建图例句柄 - 使用Rectangle而不是Patch
    legend_handles = []
    for i in range(len(legend_labels) + 1):
        # 创建带边框的彩色矩形
        rect = Rectangle((0, 0), 1, 1,
                         facecolor=colors[i],  # 填充颜色
                         edgecolor='black',  # 边框颜色
                         linewidth=1,  # 边框线宽
                         alpha=0.9)  # 透明度
        legend_handles.append(rect)

    # 添加图例到ax2
    ax2.legend(handles=legend_handles,  # 这里传入Rectangle对象列表
               labels=legend_labels,  # 对应的标签列表
               bbox_to_anchor=(1.02, 1),
               loc='upper left',
               title=f'Genus Name(Top {len(legend_labels)})',
               fontsize=8,
               title_fontsize=9,
               ncol=2,
               frameon=True,  # 确保图例外框显示
               framealpha=0.9,  # 图例外框透明度
               edgecolor='gray')  # 图例外框颜色

    # legend_labels = genus_name[:len(genus_name)//7]
    # legend_handles = [mpatches.Patch(color=colors[i], label=legend_labels[i],edgecolor='black',
    #                              linewidth=1) for i in range(len(genus_name)//7)]
    #
    # ax2.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc='upper left',
    #            title='Genus Name', fontsize=6, title_fontsize=9,ncol=2)

    # ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title='Genus Name', fontsize=8)

    total_counts = matrix_array.sum(axis=1)
    for i, total in enumerate(total_counts):
        ax2.text(i, total + 0.1, f'{int(total)}',
                 ha='center', va='bottom', fontsize=6, fontweight='bold')



    ###
    # 图表3：地理信息树状图，应当需要截取
    geo_group = csv_df['geography'].value_counts()
    if len(geo_group) > 38:
        geo_group = geo_group.head(38)

    #对na计数
    na_num = geo_group['na']
    geo_group = geo_group.drop('na')


    geo_group_index = geo_group.index.tolist()
    geo_group_value = geo_group.values
    geo_df = pd.DataFrame({'geo_name': geo_group_index, 'geo_value': geo_group_value})
    print(geo_df)

    #应当换用不同的风格颜色或者表达
    colors = get_colors(len(geo_group))

    # 设置图表3边框 - 只保留XY轴
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_linewidth(2)
    ax3.spines['bottom'].set_linewidth(2)

    #创建正确的标签列表
    tree_labels = []
    for name in geo_group_index:
        # 移除列表符号
        clean_label = name.replace("['", "").replace("']", "")

        # 检测标签是否超出区域（估算）
        # 这里我们使用一个简单的启发式方法：如果标签长度超过15个字符，则缩写
        if len(clean_label) > 11:
            if ":" in clean_label:
                # 如果有冒号，保留第一部分和最后一部分
                parts = clean_label.split(":")
                if len(parts) > 2:
                    # 多个冒号的情况
                    abbreviated_label = parts[0] + ":.." #+ parts[-1]
                else:
                    # 只有一个冒号
                    abbreviated_label = parts[0] + ":" + parts[1][:1] + ".."
            else:
                # 没有冒号，直接截断
                abbreviated_label = clean_label[:7] + ".."
        else:
            # 标签不超出，使用完整标签
            abbreviated_label = clean_label

        tree_labels.append(abbreviated_label)

    # 绘制树状图
    squarify.plot(
        sizes=geo_df['geo_value'],
        label=tree_labels,  # 使用处理后的标签
        color=colors,
        alpha=0.9,
        text_kwargs={'fontsize': 9, 'fontweight': 'bold'},
        edgecolor='black',
        linewidth=1,
        ax=ax3
    )

    # 创建带计数值的图例标签
    legend_labels = [f"{name} ({count})" for name, count in zip(geo_group_index, geo_group_value)]
    # 增加标签文字
    ax3.text(1, -2, f'{na_num} na ignored',
            horizontalalignment='left', verticalalignment='center',fontsize=8)

    # 添加图例
    patches = []
    for i, name in enumerate(geo_group_index):
        patch = plt.Rectangle((0, 0), 1, 1, fc=colors[i], edgecolor='black', linewidth=1)
        patches.append(patch)

    ax3.legend(patches, legend_labels,
               bbox_to_anchor=(1.02, 1),
               loc='upper left',
               fontsize=8,
               title=f'Geographical Location source (Top 38)')

    ax3.axis('off')  # 树状图通常不需要坐标轴
    ax3.set_title('Geographic Distribution',fontsize = 15)

    # 图表4：分类水平条形图
    list_species = csv_df['species'].value_counts()
    list_species_num = len(list_species)

    list_genus = csv_df['genus'].value_counts()
    list_genus_num = len(list_genus)

    list_family = csv_df['family'].value_counts()
    list_family_num = len(list_family)

    list_order = csv_df['order'].value_counts()
    list_order_num = len(list_order)

    list_class = csv_df['class'].value_counts()
    list_class_num = len(list_class)

    list_phylum = csv_df['phylum'].value_counts()
    list_phylum_num = len(list_phylum)

    list_kingdom = csv_df['kingdom'].value_counts()
    list_kingdom_num = len(list_kingdom)

    x_bar = [list_species_num, list_genus_num, list_family_num, list_order_num,
             list_class_num, list_phylum_num, list_kingdom_num]

    # 设置图表4边框 - 只保留XY轴
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['left'].set_linewidth(2)
    ax4.spines['bottom'].set_linewidth(2)

    bars = ax4.bar(taxonomic_levels, x_bar, width=0.5, edgecolor='black', linewidth=0.5)
    ax4.set_title('Species classification information of Reference Database',fontsize = 15)
    ax4.set_xlabel('Taxonomic Rank Level',fontsize = 12)
    ax4.set_ylabel('Number',fontsize = 12)

    for bar, value in zip(bars, x_bar):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width() / 2, height + 0.05,
                 f'{value}', ha='center', va='bottom')

    # 调整整体布局并保存
    plt.tight_layout()

    # 只保存合并后的图像
    plt.savefig(output_file,
                dpi=800, format='png', bbox_inches='tight')

