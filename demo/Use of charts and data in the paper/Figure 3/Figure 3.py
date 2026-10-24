import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.offsetbox import TextArea, HPacker, AnchoredOffsetbox

# =====================================================
# 1. 全局参数区：后续主要修改这里
# =====================================================

# ---------- 输出参数 ----------
OUTPUT_PATH = r'D:\Figure 3.png'
SAVE_DPI = 800
FIGSIZE = (12, 10)

# ---------- 字体参数 ----------
FONT_FAMILY = 'Arial'
AXIS_LABEL_SIZE = 15
TICK_LABEL_SIZE = 15
BAR_LABEL_SIZE = 13

AXIS_LABEL_WEIGHT = 'bold'

# 横坐标工具包名称加粗
XTICK_LABEL_WEIGHT = 'bold'

# 柱子顶部数字加粗
BAR_LABEL_WEIGHT = 'bold'

TEXT_COLOR = 'black'

# ---------- 子图标题参数 ----------
# 标题整体字号（调大）
TITLE_SIZE = 17

# ABCD字母加粗
PANEL_LETTER_WEIGHT = 'bold'

# 英文标题不加粗
PANEL_TEXT_WEIGHT = 'normal'

# 字母与标题之间的间距（points）
TITLE_GAP = 6

# 标题整体向右偏移：避免从y轴位置直接开始
TITLE_X_OFFSET = -0.17

# 标题整体向上偏移：拉开与子图的距离
TITLE_Y_OFFSET = 1.08

# ---------- 四张子图标题 ----------
PANEL_TITLES = {
    'A': 'Number of retrieved sequences',
    'B': 'Number of covered species',
    'C': 'Sequence retrieval efficiency',
    'D': 'Memory efficiency'
}

# ---------- 图形样式参数 ----------
BAR_WIDTH = 0.65
BAR_EDGE_COLOR = 'black'

# 普通柱子边框粗细
BAR_EDGE_LINEWIDTH = 0.8

# MFRefData柱子边框加粗参数
HIGHLIGHT_TOOL = 'MFRefData'
HIGHLIGHT_EDGE_LINEWIDTH = 1.6

SPINE_LINEWIDTH = 1.5
TICK_WIDTH = 1.2
TICK_LENGTH = 4

# ---------- 四张子图颜色参数 ----------
PANEL_COLORS = {
    'A': '#DAE3F3',
    'B': '#FAE4E3',
    'C': '#E1F2E8',
    'D': '#FFF2CC'
}

# ---------- 工具显示顺序 ----------
TOOLS = ['MFRefData', 'refdb', 'CRABS']

# ---------- 坐标轴标签 ----------
Y_LABEL_TIME_EFFICIENCY = 'Sequence Retrieval Efficiency (sequences/s)'
Y_LABEL_RAM_EFFICIENCY = 'Memory Efficiency (sequences/MB)'
Y_LABEL_SEQUENCE = 'Sequence Number'
Y_LABEL_SPECIES = 'Species Number'

# ---------- y轴刻度参数 ----------
Y_TICK_NUM = 5
Y_LIM_PAD_RATIO = 1.15

# ---------- 柱子顶部数值标签参数 ----------
BAR_LABEL_OFFSET_RATIO = 0.025

# 柱子顶部数字：整数不保留小数；非整数保留两位小数
BAR_LABEL_DECIMALS = 2

# ---------- y轴刻度小数位数 ----------
DECIMALS_TIME_EFFICIENCY = 3
DECIMALS_RAM_EFFICIENCY = 3
DECIMALS_SEQUENCE = 0
DECIMALS_SPECIES = 0

# ---------- x轴标签参数 ----------
XTICK_ROTATION = 0

# ---------- 是否打印计算结果 ----------
PRINT_CALCULATION_TABLE = True


# =====================================================
# 2. 原始数据区
# =====================================================

download_time = {
    'MFRefData': 1410.50,
    'refdb': 1285.246,
    'CRABS': 1634.3
}

processing_time = {
    'MFRefData': 2438,
    'refdb': 1.721,
    'CRABS': 520.01
}

download_ram = {
    'MFRefData': 91.1,
    'refdb': 310.892,
    'CRABS': 108.6
}

processing_ram = {
    'MFRefData': 341.47,
    'refdb': 348.753,
    'CRABS': 112.4
}

sequence_number = {
    'MFRefData': 5085,
    'refdb': 1813,
    'CRABS': 363
}

species_number = {
    'MFRefData': 725,
    'refdb': 196,
    'CRABS': 161
}


# =====================================================
# 3. 指标计算区
# =====================================================

total_time = {}
weighted_ram = {}

download_time_per_sequence = {}
total_time_per_sequence = {}
weighted_ram_per_sequence = {}

sequence_retrieval_efficiency = {}
memory_efficiency = {}

for tool in TOOLS:
    total_time[tool] = download_time[tool] + processing_time[tool]

    weighted_ram[tool] = (
        download_time[tool] * download_ram[tool]
        + processing_time[tool] * processing_ram[tool]
    ) / total_time[tool]

    # 原始指标：数值越低越好
    download_time_per_sequence[tool] = download_time[tool] / sequence_number[tool]
    total_time_per_sequence[tool] = total_time[tool] / sequence_number[tool]
    weighted_ram_per_sequence[tool] = weighted_ram[tool] / sequence_number[tool]

    # 转换后的正向效率指标：数值越高越好
    sequence_retrieval_efficiency[tool] = sequence_number[tool] / download_time[tool]
    memory_efficiency[tool] = sequence_number[tool] / weighted_ram[tool]


# =====================================================
# 4. 工具函数区
# =====================================================

def format_number(value, decimals=3):
    """用于y轴刻度：整数可不显示小数。"""
    value = float(value)

    if decimals == 0:
        return str(int(round(value)))

    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))

    return f'{value:.{decimals}f}'.rstrip('0').rstrip('.')


def format_bar_label(value, decimals=2):
    """
    用于柱子顶部数字：
    - 整数显示为整数
    - 非整数保留两位小数
    """
    value = float(value)

    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))

    return f'{value:.{decimals}f}'


def nice_ceiling(value):
    """自动生成较美观的线性y轴上限。"""
    if value <= 0:
        return 1

    magnitude = 10 ** np.floor(np.log10(value))
    multipliers = [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 10]

    for m in multipliers:
        candidate = m * magnitude
        if candidate >= value:
            return candidate

    return 10 * magnitude


def set_five_y_ticks(ax, max_value, decimals):
    """每张图固定5个y轴主刻度。"""
    upper = nice_ceiling(max_value * Y_LIM_PAD_RATIO)
    ticks = np.linspace(0, upper, Y_TICK_NUM)

    ax.set_ylim(0, upper)
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(
        FixedFormatter([format_number(t, decimals) for t in ticks])
    )


def style_axis(ax, ylabel):
    """统一坐标轴样式。"""
    ax.set_ylabel(
        ylabel,
        fontsize=AXIS_LABEL_SIZE,
        fontweight=AXIS_LABEL_WEIGHT,
        color=TEXT_COLOR
    )

    ax.tick_params(
        axis='both',
        labelsize=TICK_LABEL_SIZE,
        width=TICK_WIDTH,
        length=TICK_LENGTH,
        colors=TEXT_COLOR
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.spines['left'].set_linewidth(SPINE_LINEWIDTH)
    ax.spines['bottom'].set_linewidth(SPINE_LINEWIDTH)

    ax.spines['left'].set_color(TEXT_COLOR)
    ax.spines['bottom'].set_color(TEXT_COLOR)

    for label in ax.get_xticklabels():
        label.set_fontweight(XTICK_LABEL_WEIGHT)
        label.set_color(TEXT_COLOR)


def set_panel_title(ax, panel_id):
    """
    在每个子图顶部添加标题：
    - 只有 A./B./C./D. 加粗
    - 标题整体右移一点，不从y轴直接开始
    - 标题上移一点，和子图拉开距离
    """

    letter = TextArea(
        f'{panel_id}.',
        textprops={
            'fontfamily': FONT_FAMILY,
            'fontsize': TITLE_SIZE,
            'fontweight': PANEL_LETTER_WEIGHT,
            'color': TEXT_COLOR
        }
    )

    title_text = TextArea(
        PANEL_TITLES[panel_id],
        textprops={
            'fontfamily': FONT_FAMILY,
            'fontsize': TITLE_SIZE,
            'fontweight': PANEL_TEXT_WEIGHT,
            'color': TEXT_COLOR
        }
    )

    title_box = HPacker(
        children=[letter, title_text],
        align='baseline',
        pad=0,
        sep=TITLE_GAP
    )

    anchored_title = AnchoredOffsetbox(
        loc='lower left',
        child=title_box,
        pad=0,
        borderpad=0,
        frameon=False,
        bbox_to_anchor=(TITLE_X_OFFSET, TITLE_Y_OFFSET),
        bbox_transform=ax.transAxes
    )

    ax.add_artist(anchored_title)


def plot_single_bar(ax, values_dict, ylabel, decimals, bar_color):
    """绘制单组柱状图。"""
    values = [values_dict[tool] for tool in TOOLS]
    x = np.arange(len(TOOLS))

    bars = ax.bar(
        x,
        values,
        width=BAR_WIDTH,
        color=bar_color,
        edgecolor=BAR_EDGE_COLOR,
        linewidth=BAR_EDGE_LINEWIDTH
    )

    for bar, tool in zip(bars, TOOLS):
        if tool == HIGHLIGHT_TOOL:
            bar.set_linewidth(HIGHLIGHT_EDGE_LINEWIDTH)
            bar.set_edgecolor(BAR_EDGE_COLOR)
        else:
            bar.set_linewidth(BAR_EDGE_LINEWIDTH)
            bar.set_edgecolor(BAR_EDGE_COLOR)

    max_value = max(values)

    set_five_y_ticks(
        ax=ax,
        max_value=max_value,
        decimals=decimals
    )

    for bar, value in zip(bars, values):
        height = bar.get_height()

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + max_value * BAR_LABEL_OFFSET_RATIO,
            format_bar_label(value, BAR_LABEL_DECIMALS),
            ha='center',
            va='bottom',
            fontsize=BAR_LABEL_SIZE,
            fontweight=BAR_LABEL_WEIGHT,
            color=TEXT_COLOR
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        TOOLS,
        rotation=XTICK_ROTATION,
        color=TEXT_COLOR,
        fontweight=XTICK_LABEL_WEIGHT
    )

    style_axis(ax, ylabel)


def print_calculation_table():
    """打印计算结果，方便核对。"""
    print('\nCalculation results:')
    print('-' * 135)
    print(
        f'{"Tool":<12}'
        f'{"Download time":>16}'
        f'{"Total time":>14}'
        f'{"Weighted RAM":>16}'
        f'{"Seq number":>14}'
        f'{"Species":>10}'
        f'{"Download/Seq":>16}'
        f'{"Total/Seq":>14}'
        f'{"RAM/Seq":>12}'
        f'{"Seq/s":>12}'
        f'{"Seq/MB":>12}'
    )
    print('-' * 135)

    for tool in TOOLS:
        print(
            f'{tool:<12}'
            f'{download_time[tool]:>16.3f}'
            f'{total_time[tool]:>14.3f}'
            f'{weighted_ram[tool]:>16.3f}'
            f'{sequence_number[tool]:>14}'
            f'{species_number[tool]:>10}'
            f'{download_time_per_sequence[tool]:>16.3f}'
            f'{total_time_per_sequence[tool]:>14.3f}'
            f'{weighted_ram_per_sequence[tool]:>12.3f}'
            f'{sequence_retrieval_efficiency[tool]:>12.3f}'
            f'{memory_efficiency[tool]:>12.3f}'
        )

    print('-' * 135)


# =====================================================
# 5. 绘图主程序
# =====================================================

plt.rcParams['font.family'] = FONT_FAMILY
plt.rcParams['axes.unicode_minus'] = False

if PRINT_CALCULATION_TABLE:
    print_calculation_table()

fig, axes = plt.subplots(
    2,
    2,
    figsize=FIGSIZE,
    constrained_layout=True
)
fig.set_constrained_layout_pads(
    hspace=0.05,  # 上下两排间距
    wspace=0.05   # 左右两列间距
)
ax1, ax2, ax3, ax4 = axes.flat

# ---------- A ----------
plot_single_bar(
    ax=ax1,
    values_dict=sequence_number,
    ylabel=Y_LABEL_SEQUENCE,
    decimals=DECIMALS_SEQUENCE,
    bar_color=PANEL_COLORS['C']
)
set_panel_title(ax1, 'A')

# ---------- B ----------
plot_single_bar(
    ax=ax2,
    values_dict=species_number,
    ylabel=Y_LABEL_SPECIES,
    decimals=DECIMALS_SPECIES,
    bar_color=PANEL_COLORS['D']
)
set_panel_title(ax2, 'B')

# ---------- C ----------
plot_single_bar(
    ax=ax3,
    values_dict=sequence_retrieval_efficiency,
    ylabel=Y_LABEL_TIME_EFFICIENCY,
    decimals=DECIMALS_TIME_EFFICIENCY,
    bar_color=PANEL_COLORS['A']
)
set_panel_title(ax3, 'C')

# ---------- D ----------
plot_single_bar(
    ax=ax4,
    values_dict=memory_efficiency,
    ylabel=Y_LABEL_RAM_EFFICIENCY,
    decimals=DECIMALS_RAM_EFFICIENCY,
    bar_color=PANEL_COLORS['B']
)
set_panel_title(ax4, 'D')

plt.savefig(
    OUTPUT_PATH,
    dpi=SAVE_DPI,
    bbox_inches='tight'
)

plt.show()