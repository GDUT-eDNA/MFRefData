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

# ---------- 输入输出路径 ----------
INPUT_FILE = r"D:\Filter.csv"
OUTPUT_PATH = r'D:\Figure 4.png'  # B/C/D 组合图输出路径
SVG_OUTPUT_PATH = os.path.splitext(OUTPUT_PATH)[0] + '.svg'  # 组合图矢量输出
A_OUTPUT_PATH = os.path.splitext(OUTPUT_PATH)[0] + '_A.png'  # 原 A 子图单独输出
A_SVG_OUTPUT_PATH = os.path.splitext(A_OUTPUT_PATH)[0] + '.svg'
SAVE_DPI = 800

# ---------- 所有需要截取的地方统一前 20 ----------
TOP_N = 20

# ---------- 地理信息参数 ----------
GEO_TOP_N = TOP_N
GEO_LABEL_MAX_LEN = 8

FONT_FAMILY = 'Arial'
TITLE_SIZE = 23                  # 面板标题字号
AXIS_LABEL_SIZE = 21             # 横、纵轴名称
TICK_LABEL_SIZE = 19             # 常规刻度
SMALL_TICK_SIZE = 18             # B 图分类刻度
BAR_TEXT_SIZE = 14               # A 图堆叠柱中的数字与注释
LEGEND_SIZE = 17                 # B/D 图图例
LEGEND_TITLE_SIZE = 17           # B/D 图图例标题
TREE_TEXT_SIZE = 17              # D 图树状图内部标签

# ---------- 柱状图顶部数字 ----------
PLOT2_BAR_VALUE_SIZE = 18        # B 图柱顶总数

# ---------- 图形参数 ----------
FIGSIZE = (18, 23)               # 竖版组合图画布尺寸
FIGSIZE_A = (16, 10)             # 原 A 子图单独输出时的画布大小
SPINE_LINEWIDTH = 2
BAR_EDGE_LINEWIDTH = 0.5
BAR_ALPHA = 0.9


PLOT1_TITLE_TEMPLATE = 'Sequence Composition Across Taxonomic Ranks (Colored Top {top_n})'
PLOT1_YLABEL = 'Number of Sequence'

PLOT2_TITLE = 'Species Distribution Across Order and Genus'
PLOT2_YLABEL = 'Number of Species'

PLOT3_TITLE = 'Taxonomic Coverage of the Reference Database'
PLOT3_YLABEL = 'Number'

PLOT4_TITLE = 'Geographic Origin of Reference Sequences'

PLOT2_LEGEND_TITLE_TEMPLATE = 'Genus Name (Top {top_n})'
PLOT4_LEGEND_TITLE_TEMPLATE = 'Location Source (Top {top_n})'

# ---------- 组合图位置参数：[left, bottom, width, height] ----------
# BOTTOM_PANEL_BOUNDS 的第二项由 0.07 提高到 0.115，缩小上下 A/B 间距。
TOP_PANEL_BOUNDS = (0.09, 0.57, 0.70, 0.33)
BOTTOM_PANEL_BOUNDS = (0.09, 0.115, 0.70, 0.33)

# ---------- 内嵌 C 图参数 ----------
C_INSET_BOUNDS = (0.39, 0.28, 0.60, 0.64)  # 原 C 图下移 0.06，大小不变
C_TITLE_SIZE = 20
C_AXIS_LABEL_SIZE = 16
C_TICK_LABEL_SIZE = 15             # C 图 Y 轴刻度字体
C_X_TICK_LABEL_SIZE = 19           # C 图 X 轴分类字体，进一步放大
C_BAR_VALUE_SIZE = 14

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
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = FONT_FAMILY
plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
plt.rcParams['svg.fonttype'] = 'path'


# =====================================================
# 1. 工具函数（保留原统计与转换逻辑）
# =====================================================

def blend_with_white(hex_color, white_ratio=0.0):
    rgb = np.array(to_rgb(hex_color))
    white = np.array([1, 1, 1])
    new_rgb = rgb * (1 - white_ratio) + white * white_ratio
    return tuple(new_rgb)


def make_palette(n, base_colors, white_ratios):
    n_base = len(base_colors)
    return [
        blend_with_white(
            base_colors[i % n_base],
            white_ratios[(i // n_base) % len(white_ratios)]
        )
        for i in range(n)
    ]


def get_colors(n):
    return make_palette(n, PASTEL_BASE_COLORS, [0.00, 0.06, 0.12, 0.18, 0.24])


def get_blue_colors(n):
    return make_palette(n, BLUE_BASE_COLORS, [0.00, 0.04, 0.08, 0.12])


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
    """返回前 N 个分类的序列数、其余序列数及其余分类数。"""
    counts = csv_df[col_name].value_counts()
    top_values = counts.iloc[:matrix_len].tolist()
    top_values.extend([0] * (matrix_len - len(top_values)))
    return (
        top_values,
        int(counts.iloc[matrix_len:].sum()),
        max(0, len(counts) - matrix_len)
    )


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
    ax.set_title(
        rf'$\mathbf{{{letter}.}}$ {title}',
        fontsize=TITLE_SIZE,
        fontweight='normal',
        loc='left',
        pad=14
    )


def plot_panel_A(ax1, csv_df):
    """原 A 图，单独导出；仅删除 X 轴标题，保留分类刻度。"""
    taxonomic_levels = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    matrix_len = min(TOP_N, csv_df['species'].nunique())

    # 用统一循环获取 7 个分类层级，避免重复粘贴的统计代码。
    data = [get_top_row(csv_df, rank, matrix_len) for rank in taxonomic_levels]
    matrix_array = np.asarray([item[0] for item in data], dtype=int)
    appendix_row = np.asarray([item[1] for item in data], dtype=int)
    remained_num = [item[2] for item in data]

    colors = get_colors(matrix_len)
    bottom = np.zeros(len(taxonomic_levels))
    set_xy_spines(ax1)
    set_panel_title(ax1, 'A', PLOT1_TITLE_TEMPLATE.format(top_n=matrix_len))
    # 不再显示原来的 "Taxonomic Level" X 轴标题。
    ax1.set_ylabel(PLOT1_YLABEL, fontsize=AXIS_LABEL_SIZE)
    ax1.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)

    for i in range(matrix_len):
        values = matrix_array[:, i]
        ax1.bar(
            taxonomic_levels, values, bottom=bottom,
            color=colors[i], alpha=BAR_ALPHA,
            edgecolor='black', linewidth=0.3
        )
        bottom += values

    totals = bottom.copy()
    for i, total in enumerate(totals):
        cumulative = 0
        for value in matrix_array[i]:
            if total > 0 and value > 0 and value / total > 0.009:
                ax1.text(
                    i, cumulative + value / 2, f'{int(value)}',
                    ha='center', va='center', fontsize=BAR_TEXT_SIZE,
                    fontweight='bold', color='black'
                )
            cumulative += value

    ax1.bar(
        taxonomic_levels, appendix_row, bottom=totals,
        color='whitesmoke', alpha=0.8,
        edgecolor='black', linewidth=0.3
    )

    for i, remaining in enumerate(appendix_row):
        if remaining == 0:
            continue
        ax1.text(
            i, totals[i] + remaining / 2,
            f'remaining {int(remaining)} seqs\nbelonged to {remained_num[i]} {taxonomic_levels[i]}',
            ha='center', va='center', fontsize=BAR_TEXT_SIZE,
            fontweight='bold', color='maroon'
        )


def plot_panel_B_with_inset(ax2, ax4, csv_df):
    order_name_label = csv_df['order'].value_counts().index.tolist()
    genus_name_value = csv_df['genus'].value_counts()
    genus_name = genus_name_value.index.tolist()

    present = pd.crosstab(csv_df['order'], csv_df['genus']).reindex(
        index=order_name_label, columns=genus_name, fill_value=0
    ).gt(0)

    species_per_genus = csv_df.groupby('genus')['species'].nunique()
    matrix_array = present.mul(species_per_genus.reindex(genus_name), axis=1).to_numpy(dtype=int)
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

    set_panel_title(ax2, 'A', PLOT2_TITLE)
    ax2.set_ylabel(PLOT2_YLABEL, fontsize=AXIS_LABEL_SIZE)
    ax2.tick_params(axis='y', labelsize=TICK_LABEL_SIZE)
    ax2.set_xticks(range(len(order_name_label)))

    order_display_labels = [name[2:] if name.startswith('o_') else name for name in order_name_label]
    ax2.set_xticklabels(
        order_display_labels,
        rotation=90,
        ha='center',
        va='top',
        fontsize=SMALL_TICK_SIZE,
        fontweight='normal'
    )
    ax2.tick_params(axis='x', pad=5, length=4)

    legend_labels = genus_name[:min(TOP_N, len(genus_name))]
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=colors[i], edgecolor='black', linewidth=1, alpha=BAR_ALPHA)
        for i in range(len(legend_labels))
    ]
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

    # -------- inset C：直接覆盖在原 B 图右上角空白区域 --------
    ax4.set_facecolor('white')
    ax4.patch.set_alpha(1.0)
    set_xy_spines(ax4)
    ax4.set_zorder(10)

    taxonomic_levels = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    x_bar = [
        csv_df['species'].nunique(),
        csv_df['genus'].nunique(),
        csv_df['family'].nunique(),
        csv_df['order'].nunique(),
        csv_df['class'].nunique(),
        csv_df['phylum'].nunique(),
        csv_df['kingdom'].nunique()
    ]

    bar_colors = get_blue_colors(len(x_bar))
    inset_x = np.arange(len(taxonomic_levels))
    bars = ax4.bar(
        inset_x, x_bar, width=0.5,
        color=bar_colors,
        edgecolor='black',
        linewidth=BAR_EDGE_LINEWIDTH
    )

    ax4.set_title(PLOT3_TITLE, fontsize=C_TITLE_SIZE, fontweight='normal', loc='left', pad=8)
    ax4.set_ylabel(PLOT3_YLABEL, fontsize=C_AXIS_LABEL_SIZE)
    ax4.tick_params(axis='y', labelsize=C_TICK_LABEL_SIZE)
    ax4.tick_params(axis='x', labelsize=C_X_TICK_LABEL_SIZE)
    ax4.set_xticks(inset_x)
    ax4.set_xticklabels(taxonomic_levels, rotation=90, ha='center')
    ax4.set_xlim(-0.6, len(taxonomic_levels) - 0.4)

    y_max = max(x_bar) if x_bar else 1
    ax4.set_ylim(0, y_max * 1.18)

    for bar_obj, value in zip(bars, x_bar):
        height = bar_obj.get_height()
        ax4.text(
            bar_obj.get_x() + bar_obj.get_width() / 2,
            height + max(y_max * 0.015, 0.05), f'{value}',
            ha='center', va='bottom',
            fontsize=C_BAR_VALUE_SIZE,
            fontweight='bold', color='black'
        )


def plot_panel_D(ax3, csv_df):
    if 'geography' not in csv_df.columns:
        print('Warning: geography column not found.')
        geo_group = pd.Series(dtype=int)
    else:
        geo_group_all = csv_df['geography'].map(normalize_geo_to_country).value_counts()
        geo_group = geo_group_all.drop(labels=['na'], errors='ignore').head(GEO_TOP_N)

    geo_group_index = geo_group.index.tolist()
    geo_group_value = geo_group.values
    geo_df = pd.DataFrame({'geo_name': geo_group_index, 'geo_value': geo_group_value})
    colors = get_colors(len(geo_group))

    if len(geo_df) > 0:
        geo_full_labels = [clean_display_label(name) for name in geo_group_index]
        legend_labels = [f'{name} ({count})' for name, count in zip(geo_full_labels, geo_group_value)]
        tree_labels = [geo_label_filter(name, GEO_LABEL_MAX_LEN) for name in geo_full_labels]

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

        patches = [plt.Rectangle((0, 0), 1, 1, fc=colors[i], edgecolor='black', linewidth=1)
                   for i in range(len(geo_group_index))]
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

    set_panel_title(ax3, 'B', PLOT4_TITLE)
    set_treemap_axis_style(ax3)


def save_png_svg(fig, png_path, svg_path):
    """同一版布局保存为高分辨率 PNG 与可缩放 SVG。"""
    for filepath in (png_path, svg_path):
        os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
        fig.savefig(filepath, dpi=SAVE_DPI, bbox_inches='tight')
        print(f'已保存：{os.path.abspath(filepath)}')


def bar(file):
    """读取一次输入数据，生成单独 A 图和组合 A/B 图。"""
    if not os.path.isfile(file):
        print(f'File not found: {file}')
        return

    try:
        csv_df = pd.read_csv(file, encoding='utf-8-sig')
    except Exception as e:
        print(f'读取 CSV 失败：{e}')
        return

    required_columns = {'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom'}
    missing_columns = required_columns.difference(csv_df.columns)
    if missing_columns:
        print('缺少必要字段：' + ', '.join(sorted(missing_columns)))
        return

    # A：原 A 子图单独输出，不带 "Taxonomic Level"。
    fig_a, ax_a = plt.subplots(figsize=FIGSIZE_A, constrained_layout=True)
    try:
        plot_panel_A(ax_a, csv_df)
        save_png_svg(fig_a, A_OUTPUT_PATH, A_SVG_OUTPUT_PATH)
    finally:
        plt.close(fig_a)

    # 组合图上方为 A（原 B + 内嵌 C），下方为 B（原 D）。
    fig = plt.figure(figsize=FIGSIZE)
    try:
        ax_top = fig.add_axes(TOP_PANEL_BOUNDS)
        ax_bottom = fig.add_axes(BOTTOM_PANEL_BOUNDS)
        ax_inset = ax_top.inset_axes(C_INSET_BOUNDS)

        plot_panel_B_with_inset(ax_top, ax_inset, csv_df)
        plot_panel_D(ax_bottom, csv_df)
        save_png_svg(fig, OUTPUT_PATH, SVG_OUTPUT_PATH)
    finally:
        plt.close(fig)


if __name__ == '__main__':
    bar(INPUT_FILE)
