import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import squarify
import re
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb


# =====================================================
# 0. 全局参数区：后续主要修改这里
# =====================================================

# ---------- 输入输出路径（保留原路径，使用前请确认） ----------
INPUT_FILE = r'D:\Filter.csv'
OUTPUT_PATH = r'D:\Figure 4.png'
SVG_OUTPUT_PATH = os.path.splitext(OUTPUT_PATH)[0] + '.svg'  # 与 PNG 同目录、同名
SAVE_DPI = 800

# ---------- 所有需要截取的地方统一前 20 ----------
TOP_N = 20

# ---------- 地理信息参数 ----------
GEO_TOP_N = TOP_N
GEO_LABEL_MAX_LEN = 8

# ---------- 字体参数：统一、易读，同时保留各元素间的字号层次 ----------
FONT_FAMILY = 'Arial'
TITLE_SIZE = 22                  # 四个子图标题：仅 A./B./C./D. 加粗
AXIS_LABEL_SIZE = 20             # 横、纵轴名称：略微放大
TICK_LABEL_SIZE = 18             # 常规刻度
SMALL_TICK_SIZE = 17             # B 图分类较密，避免字体过大导致重叠
BAR_TEXT_SIZE = 13               # A 图堆叠柱中的数字与注释
LEGEND_SIZE = 16                 # B/D 图图例
LEGEND_TITLE_SIZE = 16           # B/D 图图例标题
TREE_TEXT_SIZE = 16              # D 图树状图内部标签
NA_TEXT_SIZE = 15                # D 图 na 数量注释

# ---------- 柱状图顶部数字 ----------
PLOT2_BAR_VALUE_SIZE = 17        # B 图柱顶总数
PLOT3_BAR_VALUE_SIZE = 17        # C 图柱顶总数

# ---------- 图形参数 ----------
FIGSIZE = (27, 19)
SPINE_LINEWIDTH = 2
BAR_EDGE_LINEWIDTH = 0.5
BAR_ALPHA = 0.9

# =====================================================
# 0.1 四个输出位置的标题与坐标轴名称
# 左上 A: ax1 | 右上 B: ax2 | 左下 C: ax4 | 右下 D: ax3
# =====================================================

PLOT1_TITLE_TEMPLATE = 'Sequence Composition Across Taxonomic Ranks (Colored Top {top_n})'
PLOT1_XLABEL = 'Taxonomic Level'
PLOT1_YLABEL = 'Number of Sequence'

PLOT2_TITLE = 'Species Distribution Across Order and Genus'
PLOT2_XLABEL = 'Order'
PLOT2_YLABEL = 'Number of Species'

PLOT3_TITLE = 'Taxonomic Coverage of the Reference Database'
PLOT3_XLABEL = 'Taxonomic Rank Level'
PLOT3_YLABEL = 'Number'

PLOT4_TITLE = 'Geographic Origin of Reference Sequences'
PLOT4_XLABEL = ''
PLOT4_YLABEL = ''

PLOT2_LEGEND_TITLE_TEMPLATE = 'Genus Name (Top {top_n})'
PLOT4_LEGEND_TITLE_TEMPLATE = 'Location Source (Top {top_n})'

# ---------- 原始颜色定义：不改变颜色映射 ----------
PASTEL_BASE_COLORS = [
    '#ADBBD9', '#E6A19A', '#A3D4BA', '#F2D7A2',
    '#B7DDE0', '#DFAEAA', '#C2D8EE', '#C6DDB0',
    '#D8B7DD', '#E7B8B3', '#D5C6A8', '#D9D9D9'
]

BLUE_BASE_COLORS = [
    '#3F7FB8', '#4F8FC9', '#5D9AD1', '#6BA6D9', '#7AB1DF',
    '#3472AA', '#2F6699', '#5A8FBE', '#6F9FCA', '#83AED4'
]

plt.rcParams['font.family'] = FONT_FAMILY
plt.rcParams['font.sans-serif'] = [FONT_FAMILY]
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.titlepad'] = 14
# 数学文本与正文尽量匹配 Arial，使仅加粗的面板字母不显突兀
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = FONT_FAMILY
plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
# 将 SVG 字体保存为真正的矢量轮廓，避免其他电脑缺少 Arial 时字体被替换
plt.rcParams['svg.fonttype'] = 'path'


# =====================================================
# 1. 工具函数（保留原统计与转换逻辑）
# =====================================================

def blend_with_white(hex_color, white_ratio=0.0):
    rgb = np.array(to_rgb(hex_color))
    white = np.array([1, 1, 1])
    new_rgb = rgb * (1 - white_ratio) + white * white_ratio
    return tuple(new_rgb)


def get_colors(n):
    colors = []
    for i in range(n):
        base_color = PASTEL_BASE_COLORS[i % len(PASTEL_BASE_COLORS)]
        cycle = i // len(PASTEL_BASE_COLORS)
        white_ratio_list = [0.00, 0.06, 0.12, 0.18, 0.24]
        white_ratio = white_ratio_list[cycle % len(white_ratio_list)]
        colors.append(blend_with_white(base_color, white_ratio))
    return colors


def get_blue_colors(n):
    colors = []
    for i in range(n):
        base_color = BLUE_BASE_COLORS[i % len(BLUE_BASE_COLORS)]
        cycle = i // len(BLUE_BASE_COLORS)
        white_ratio_list = [0.00, 0.04, 0.08, 0.12]
        white_ratio = white_ratio_list[cycle % len(white_ratio_list)]
        colors.append(blend_with_white(base_color, white_ratio))
    return colors


def clean_display_label(label):
    if pd.isna(label):
        return 'na'
    label = str(label).strip()
    label = label.replace('[', '').replace(']', '')
    label = label.replace("'", '').replace('"', '')
    label = re.sub(r'\s+', ' ', label)
    return label.strip()


def truncate_label(label, max_len=15):
    label = clean_display_label(label)
    if len(label) <= max_len:
        return label
    if max_len <= 2:
        return label[:max_len]
    return label[:max_len - 2] + '..'


def normalize_geo_to_country(value):
    if pd.isna(value):
        return 'na'

    text = clean_display_label(value)
    if text == '':
        return 'na'
    low = text.lower().strip()
    if low in ['na', 'nan', 'none', 'null', 'missing', 'not available', 'not applicable']:
        return 'na'

    text = re.split(r'\s*;\s*|\s*\|\s*', text)[0].strip()
    if text == '':
        return 'na'
    if ':' in text:
        country = text.split(':', 1)[0].strip()
    else:
        country = text.split(',', 1)[0].strip()
    if country == '':
        return 'na'

    country = clean_display_label(country)
    country_low = country.lower().strip()

    china_aliases = {
        'china', 'p. r. china', 'pr china', 'p.r. china',
        'people’s republic of china', "people's republic of china",
        'peoples republic of china', 'mainland china', 'taiwan',
        'taiwan province', 'taiwan province of china', 'chinese taipei',
        'republic of china', 'roc', 'hong kong', 'hongkong',
        'hong kong sar', 'hong kong sar china',
        'hong kong special administrative region', 'macau', 'macao',
        'macau sar', 'macao sar', 'macao sar china', 'macau sar china',
        'macau special administrative region',
        'macao special administrative region', 'tibet', 'xizang', 'xinjiang'
    }
    if country_low in china_aliases:
        return 'China'

    if (
        'china' in country_low or 'taiwan' in country_low or
        'hong kong' in country_low or 'hongkong' in country_low or
        'macau' in country_low or 'macao' in country_low or
        'chinese taipei' in country_low
    ):
        return 'China'

    country_standard = {
        'U.S.A.': 'USA', 'US': 'USA',
        'United States': 'USA', 'United States of America': 'USA'
    }
    return country_standard.get(country, country)


def geo_label_filter(label, max_len=15):
    return truncate_label(label, max_len=max_len)


def get_top_row(csv_df, col_name, matrix_len):
    row = []
    count = csv_df[col_name].value_counts()
    if len(count) >= matrix_len:
        row.extend(count.head(matrix_len).tolist())
        appendix_value = count.values[matrix_len:].sum()
        remained_number = len(count) - matrix_len
    else:
        row.extend(count.values.tolist())
        while len(row) < matrix_len:
            row.append(0)
        appendix_value = 0
        remained_number = 0
    return row, appendix_value, remained_number


def set_xy_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(SPINE_LINEWIDTH)
    ax.spines['bottom'].set_linewidth(SPINE_LINEWIDTH)


def set_treemap_axis_style(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def set_panel_title(ax, letter, title):
    """标题与子图左上角对齐，仅加粗面板编号，不加粗英文标题。"""
    ax.set_title(
        rf'$\mathbf{{{letter}.}}$ {title}',
        fontsize=TITLE_SIZE,
        fontweight='normal',
        loc='left',              # 从子图坐标区域的左边界开始，不居中
        pad=14                   # 标题和子图上边界留出适当空隙
    )


# =====================================================
# 2. 主绘图函数
# =====================================================

def bar(file):
    print('plotting')

    if os.path.exists(file):
        try:
            csv_df = pd.read_csv(file, encoding='utf-8-sig')
        except Exception as e:
            print(e)
            return
    else:
        print(f'File not found: {file}')
        return

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE, constrained_layout=True)
    ax1, ax2, ax4, ax3 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # =====================================================
    # A. 左上 ax1：分类组成堆叠条形图
    # =====================================================
    taxonomic_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomic_levels.reverse()

    appendix_row = [0] * 7
    remained_num = [0] * 7

    species_count = csv_df['species'].value_counts()
    matrix_len = min(TOP_N, len(species_count))

    species_row, appendix_row[0], remained_num[0] = get_top_row(csv_df, 'species', matrix_len)
    genus_row, appendix_row[1], remained_num[1] = get_top_row(csv_df, 'genus', matrix_len)
    family_row, appendix_row[2], remained_num[2] = get_top_row(csv_df, 'family', matrix_len)
    order_row, appendix_row[3], remained_num[3] = get_top_row(csv_df, 'order', matrix_len)
    class_row, appendix_row[4], remained_num[4] = get_top_row(csv_df, 'class', matrix_len)
    phylum_row, appendix_row[5], remained_num[5] = get_top_row(csv_df, 'phylum', matrix_len)
    kingdom_row, appendix_row[6], remained_num[6] = get_top_row(csv_df, 'kingdom', matrix_len)

    matrix_df = [
        species_row, genus_row, family_row, order_row,
        class_row, phylum_row, kingdom_row
    ]

    colors = get_colors(matrix_len)
    bottom = np.zeros(len(taxonomic_levels))
    matrix_array = np.array(matrix_df)
    set_xy_spines(ax1)

    # 左上角对齐，只有 A. 加粗
    set_panel_title(ax1, 'A', PLOT1_TITLE_TEMPLATE.format(top_n=matrix_len))
    ax1.set_xlabel(PLOT1_XLABEL, fontsize=AXIS_LABEL_SIZE)
    ax1.set_ylabel(PLOT1_YLABEL, fontsize=AXIS_LABEL_SIZE)
    ax1.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)

    bar_heights = np.zeros((len(taxonomic_levels), matrix_len))
    bottom_last = [[]]

    for i in range(matrix_len):
        row = matrix_array[:, i]
        bars = ax1.bar(
            taxonomic_levels, row, bottom=bottom,
            color=colors[i], alpha=BAR_ALPHA,
            edgecolor='black', linewidth=0.3
        )
        for j, bar_obj in enumerate(bars):
            bar_heights[j, i] = row[j]
        bottom += row
        bottom_last = bottom

    for i in range(len(taxonomic_levels)):
        total_height = bottom_last[i]
        cumulative_height = 0
        for j in range(matrix_len):
            height = bar_heights[i, j]
            if total_height > 0 and height > 0 and height / total_height > 0.009:
                label_y = cumulative_height + height / 2
                ax1.text(
                    i, label_y, f'{int(height)}',
                    ha='center', va='center', fontsize=BAR_TEXT_SIZE,
                    fontweight='bold', color='black'
                )
            cumulative_height += height

    ax1.bar(
        taxonomic_levels, appendix_row, bottom=bottom_last,
        color='whitesmoke', alpha=0.8,
        edgecolor='black', linewidth=0.3
    )

    for i in range(len(taxonomic_levels)):
        if appendix_row[i] == 0:
            continue
        label_y = bottom_last[i] + appendix_row[i] / 2
        ax1.text(
            i, label_y,
            f'remaining {int(appendix_row[i])} seqs\nbelonged to {remained_num[i]} {taxonomic_levels[i]}',
            ha='center', va='center', fontsize=BAR_TEXT_SIZE,
            fontweight='bold', color='maroon'
        )

    # =====================================================
    # B. 右上 ax2：目-属堆叠条形图
    # =====================================================
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
            if (csv_df[csv_df['genus'] == single_genus]['order'] == single_order).any():
                single_row.append(csv_df[csv_df['genus'] == single_genus]['species'].nunique())
            else:
                single_row.append(0)
        order_matrix.append(single_row)

    matrix_array = np.array(order_matrix).T
    colors = get_colors(len(genus_name))
    bottom = np.zeros(len(order_name_label))
    set_xy_spines(ax2)

    for i, genus in enumerate(genus_name):
        counts = matrix_array[:, i]
        ax2.bar(
            order_name_label, counts, bottom=bottom,
            label=genus, color=colors[i], alpha=BAR_ALPHA,
            edgecolor='black', linewidth=BAR_EDGE_LINEWIDTH
        )
        bottom += counts

    # 左上角对齐，只有 B. 加粗
    set_panel_title(ax2, 'B', PLOT2_TITLE)
    ax2.set_xlabel(PLOT2_XLABEL, fontsize=AXIS_LABEL_SIZE)
    ax2.set_ylabel(PLOT2_YLABEL, fontsize=AXIS_LABEL_SIZE)
    ax2.tick_params(axis='y', labelsize=TICK_LABEL_SIZE)
    ax2.set_xticks(range(len(order_name_label)))

    # 修改2：只改变右上横坐标显示大小，不改变分类数据或顺序。
    # 仍仅移除真实存在的 o_ 前缀，保持垂直排列避免长名称重叠。
    order_display_labels = [
        name[2:] if name.startswith('o_') else name
        for name in order_name_label
    ]
    ax2.set_xticklabels(
        order_display_labels,
        rotation=90,
        ha='center',
        va='top',
        fontsize=SMALL_TICK_SIZE,
        fontweight='normal'      # 由 medium 调为 normal，字重与其他刻度统一
    )
    ax2.tick_params(axis='x', pad=5, length=4)

    legend_labels = genus_name[:min(TOP_N, len(genus_name))]
    legend_handles = []
    for i in range(len(legend_labels)):
        legend_handles.append(Rectangle(
            (0, 0), 1, 1, facecolor=colors[i],
            edgecolor='black', linewidth=1, alpha=BAR_ALPHA
        ))

    ax2.legend(
        handles=legend_handles,
        labels=legend_labels,
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        title=PLOT2_LEGEND_TITLE_TEMPLATE.format(top_n=len(legend_labels)),
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
        ncol=1,
        frameon=True,
        framealpha=0.9,
        edgecolor='gray'
    )

    total_counts = matrix_array.sum(axis=1)
    ax2.set_ylim(bottom=0, top=max(1, float(max(total_counts, default=0)) * 1.18))

    for i, total in enumerate(total_counts):
        ax2.text(
            i, total + 0.1, f'{int(total)}',
            ha='center', va='bottom',
            fontsize=PLOT2_BAR_VALUE_SIZE,
            fontweight='bold', color='black'
        )

    # =====================================================
    # D. 右下 ax3：地理信息树状图
    # =====================================================
    if 'geography' not in csv_df.columns:
        print('Warning: geography column not found.')
        geo_group = pd.Series(dtype=int)
        na_num = 0
    else:
        csv_df['geo_country'] = csv_df['geography'].apply(normalize_geo_to_country)
        geo_group_all = csv_df['geo_country'].value_counts()
        na_num = int(geo_group_all.get('na', 0))
        geo_group = geo_group_all.drop(labels=['na'], errors='ignore')
        geo_group = geo_group.head(GEO_TOP_N)

    geo_group_index = geo_group.index.tolist()
    geo_group_value = geo_group.values
    geo_df = pd.DataFrame({'geo_name': geo_group_index, 'geo_value': geo_group_value})
    print(geo_df)
    colors = get_colors(len(geo_group))

    if len(geo_df) > 0:
        geo_full_labels = [clean_display_label(name) for name in geo_group_index]
        legend_labels = [
            f'{name} ({count})'
            for name, count in zip(geo_full_labels, geo_group_value)
        ]
        tree_labels = [
            geo_label_filter(name, GEO_LABEL_MAX_LEN)
            for name in geo_full_labels
        ]

        print('Full geographic labels used in legend:')
        print(legend_labels)
        print('Geographic labels used in treemap:')
        print(tree_labels)

        squarify.plot(
            sizes=geo_df['geo_value'],
            label=tree_labels,
            color=colors,
            alpha=BAR_ALPHA,
            text_kwargs={
                'fontsize': TREE_TEXT_SIZE,
                'fontweight': 'bold',
                'color': 'black',
                'fontfamily': FONT_FAMILY
            },
            edgecolor='black',
            linewidth=1,
            ax=ax3
        )

        ax3.text(
            1, -2, f'{na_num} na ignored',
            horizontalalignment='left',
            verticalalignment='center',
            fontsize=NA_TEXT_SIZE,
            color='black',
            fontfamily=FONT_FAMILY
        )

        patches = []
        for i, name in enumerate(geo_group_index):
            patches.append(plt.Rectangle(
                (0, 0), 1, 1, fc=colors[i],
                edgecolor='black', linewidth=1
            ))

        ax3.legend(
            patches,
            legend_labels,
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            fontsize=LEGEND_SIZE,
            title=PLOT4_LEGEND_TITLE_TEMPLATE.format(top_n=GEO_TOP_N),
            title_fontsize=LEGEND_TITLE_SIZE,
            frameon=True,
            framealpha=0.9,
            edgecolor='gray'
        )
    else:
        ax3.text(
            0.5, 0.5, 'No geographical information',
            ha='center', va='center', fontsize=TITLE_SIZE,
            fontweight='bold', color='black', fontfamily=FONT_FAMILY
        )

    # 右下角对应 D.；只有 D. 加粗
    set_panel_title(ax3, 'D', PLOT4_TITLE)
    ax3.set_xlabel(PLOT4_XLABEL, fontsize=AXIS_LABEL_SIZE)
    ax3.set_ylabel(PLOT4_YLABEL, fontsize=AXIS_LABEL_SIZE)
    set_treemap_axis_style(ax3)

    # =====================================================
    # C. 左下 ax4：分类水平条形图（原图实际为竖向柱状图）
    # =====================================================
    list_species_num = csv_df['species'].nunique()
    list_genus_num = csv_df['genus'].nunique()
    list_family_num = csv_df['family'].nunique()
    list_order_num = csv_df['order'].nunique()
    list_class_num = csv_df['class'].nunique()
    list_phylum_num = csv_df['phylum'].nunique()
    list_kingdom_num = csv_df['kingdom'].nunique()

    x_bar = [
        list_species_num, list_genus_num, list_family_num,
        list_order_num, list_class_num, list_phylum_num,
        list_kingdom_num
    ]

    set_xy_spines(ax4)
    bar_colors = get_blue_colors(len(x_bar))
    bars = ax4.bar(
        taxonomic_levels, x_bar, width=0.5,
        color=bar_colors,
        edgecolor='black',
        linewidth=BAR_EDGE_LINEWIDTH
    )

    # 左下角对齐，只有 C. 加粗
    set_panel_title(ax4, 'C', PLOT3_TITLE)
    ax4.set_xlabel(PLOT3_XLABEL, fontsize=AXIS_LABEL_SIZE)
    ax4.set_ylabel(PLOT3_YLABEL, fontsize=AXIS_LABEL_SIZE)
    ax4.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)

    for bar_obj, value in zip(bars, x_bar):
        height = bar_obj.get_height()
        ax4.text(
            bar_obj.get_x() + bar_obj.get_width() / 2,
            height + 0.05, f'{value}',
            ha='center', va='bottom',
            fontsize=PLOT3_BAR_VALUE_SIZE,
            fontweight='bold', color='black'
        )

    # =====================================================
    # 保存：保留原分辨率和输出设置
    # =====================================================
    # 使用同一个 Figure 连续导出：不重绘、不重新计算数据、不更改布局
    plt.savefig(
        OUTPUT_PATH,
        dpi=SAVE_DPI,
        format='png',
        bbox_inches='tight'
    )
    # PNG 保存完成后冻结约束布局，避免 SVG 后端重新调整四个子图的位置。
    # 此操作不会更改已经保存的 PNG，也不会改变任何统计数据或图形元素。
    fig.set_constrained_layout(False)
    plt.savefig(
        SVG_OUTPUT_PATH,
        format='svg',
        bbox_inches='tight'
    )
    print(f'已保存 PNG：{os.path.abspath(OUTPUT_PATH)}')
    print(f'已保存 SVG：{os.path.abspath(SVG_OUTPUT_PATH)}')
    print('已应用：仅 A–D 面板编号加粗、所有标题左对齐，并适当放大较小的图内字体。')


if __name__ == '__main__':
    bar(file=INPUT_FILE)
