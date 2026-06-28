from Bio.Seq import Seq
from Bio import SeqIO
import os
import pandas as pd
import numpy as np

degenerate_base_set = {
    'A': {'A'}, 'T': {'T'}, 'C': {'C'}, 'G': {'G'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
}


# ============================================================
# 1. 快速反向互补
#    用 Python 字符串 translate，比每次调用 Seq(...).reverse_complement() 更快
# ============================================================

_RC_TABLE = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn"
)


def reverse_complement_fast(seq):
    return seq.translate(_RC_TABLE)[::-1].upper()


# ============================================================
# 2. primer 预处理
#    F/R 仍然按 5' -> 3' 输入
#    R 在这里统一转换为 reverse complement
# ============================================================

def prepare_primers(F, R):
    str_f = F.upper()
    str_r = reverse_complement_fast(R)

    return {
        "F": str_f,
        "R_rc": str_r,
        "len_f": len(str_f),
        "len_r": len(str_r)
    }


def calculate_mismatch_optimized(primer, segment, degenerate_base_set, max_mismatch):
    if len(primer) != len(segment):
        return len(primer)

    mismatch = 0
    for i in range(len(primer)):
        base_primer = primer[i]
        base_segment = segment[i]

        allowed_bases = degenerate_base_set.get(base_primer, {base_primer})

        if base_segment not in allowed_bases:
            mismatch += 1
            if mismatch > max_mismatch:  # 提前终止
                return mismatch

    return mismatch


# ============================================================
# 3. 保留原函数名，但支持两种调用方式
#
#    原调用方式仍然有效：
#    remove_primers_optimized(seq, F, R, max_mismatch=mis)
#
#    新的快速调用方式：
#    primer_info = prepare_primers(F, R)
#    remove_primers_optimized(seq, primer_info, max_mismatch=mis)
# ============================================================

def remove_primers_optimized(seq, F_or_primer_info, R=None, max_mismatch=0):
    str_seq = seq.upper()

    if isinstance(F_or_primer_info, dict):
        primer_info = F_or_primer_info
    else:
        primer_info = prepare_primers(F_or_primer_info, R)

    str_f = primer_info["F"]
    str_r = primer_info["R_rc"]
    len_f = primer_info["len_f"]
    len_r = primer_info["len_r"]

    seq_len = len(str_seq)

    if seq_len < len_f or seq_len < len_r:
        return 0

    # 查找前向引物位置
    forward_positions = []
    for i in range(0, seq_len - len_f + 1, 1):
        segment = str_seq[i:i + len_f]
        mismatch = calculate_mismatch_optimized(str_f, segment, degenerate_base_set, max_mismatch)

        if mismatch <= max_mismatch:
            forward_positions.append((i, i + len_f))
            break

    if not forward_positions:
        return 0

    # 查找反向引物位置
    reverse_positions = []
    min_start = forward_positions[0][1]

    for j in range(min_start, seq_len - len_r + 1, 1):
        segment = str_seq[j:j + len_r]
        mismatch = calculate_mismatch_optimized(str_r, segment, degenerate_base_set, max_mismatch)

        if mismatch <= max_mismatch:
            reverse_positions.append((j, j + len_r))
            break

    if not reverse_positions:
        return 0

    # 返回扩增产物
    f_start, f_end = forward_positions[0]
    r_start, r_end = reverse_positions[0]

    if r_start > f_end:
        return str_seq[f_end:r_start]

    return 0


# ============================================================
# 4. 新增：双方向搜索
#    先查原始方向；
#    只有原始方向失败时，才对整条 reference 做 reverse complement 后再查一次。
#
#    这一步解决：
#    NCBI 中部分序列反向存放时，原代码会误判为不可扩增的问题。
# ============================================================

def remove_primers_both_orientations(seq, primer_info, max_mismatch=0):
    # 1. 原始方向搜索
    new_seq_record = remove_primers_optimized(
        seq,
        primer_info,
        max_mismatch=max_mismatch
    )

    if new_seq_record != 0:
        return new_seq_record

    # 2. 原始方向失败后，再检查整条序列 reverse complement
    seq_rc = reverse_complement_fast(seq)

    new_seq_record_rc = remove_primers_optimized(
        seq_rc,
        primer_info,
        max_mismatch=max_mismatch
    )

    if new_seq_record_rc != 0:
        return new_seq_record_rc

    return 0


def isPCR_optimized(input_path, F, R, mis):
    # 为了确保通用，只要是数字文件＋fasta格式都可以去引物
    if not os.path.exists(input_path):
        return print('directory is not exist')
    if not os.path.exists(input_path + '/sequences'):
        return print('sequence files is not exist')
    if not os.path.exists(input_path + '/sequence_info.csv'):
        return print('csv file is not exist')

    # 新增：primer 只预处理一次，减少每条序列重复计算
    # 注意：F/R 均按 5' -> 3' 输入，R 不要提前反向互补
    primer_info = prepare_primers(F, R)

    # 创建输出目录
    output_dir = os.path.dirname(input_path)
    output_path = os.path.join(output_dir, 'ISPCR')
    os.makedirs(output_path, exist_ok=True)

    sequence_output_path = os.path.join(output_path, 'ISPCR_sequence')
    os.makedirs(sequence_output_path, exist_ok=True)

    files = os.listdir(input_path + '/sequences')

    # 原来是 list；这里改为 set，加快后面判断 GI 是否失败
    # 输出结果不变
    file_seqs_needmove = set()

    seq_info = []

    # 修改点1：添加一个集合来记录实际遇到的GI
    encountered_gis = set()

    # 处理每个物种目录
    for file in files:
        file_parts = file.split('_')

        if len(file_parts) == 2:
            file_part = file_parts[1]
        elif len(file_parts) == 3:
            file_part = file_parts[2]
        elif len(file_parts) > 3:
            print(file)
            continue
        else:
            print(f'{file} is not taxid file, already ignored')
            continue

        if not file_part.isdigit():
            print(f'{file} is not taxid file, already ignored')
            continue

        file_seqs = []
        files_sub = os.listdir(input_path + '/sequences/' + file)
        new_combined_species_seqs = f'{sequence_output_path}/{file}_combined.fasta'

        # 处理每个FASTA文件
        for files_single in files_sub:
            if not files_single.endswith('.fasta'):
                print(f'{files_single} is not .fasta, already ignored')
                continue

            file_path = input_path + '/sequences/' + file + '/' + files_single

            for seq_record in SeqIO.parse(file_path, 'fasta'):
                seq_length = len(seq_record.seq)

                # 核心修改：
                # 原来：只检查原始方向
                # 现在：先检查原始方向，失败后再检查 reverse-complement 方向
                new_seq_record = remove_primers_both_orientations(
                    str(seq_record.seq),
                    primer_info,
                    max_mismatch=mis
                )

                gi = files_single.split('.')[0]

                # 修改点2：记录遇到的GI
                encountered_gis.add(gi)

                if new_seq_record != 0:
                    seq_record.seq = Seq(new_seq_record)
                    seq_length_after_remove = len(seq_record.seq)
                    file_seqs.append(seq_record)
                    print(f'\t物种taxid{file}的序列{files_single}已成功扩增')
                else:
                    print(f'\t物种taxid{file}的序列{files_single}不可被扩增')
                    # 依然保存，保持你原有输出逻辑不变
                    seq_length_after_remove = 0
                    file_seqs.append(seq_record)
                    file_seqs_needmove.add(files_single.split('.')[0])

                seq_info.append({
                    "GI": files_single.split('.')[0],
                    "raw_seq_length": seq_length,
                    "amplicon_length": seq_length_after_remove
                })

        # 保存合并的序列文件
        if file_seqs:
            SeqIO.write(file_seqs, new_combined_species_seqs, 'fasta')
            print(f'物种{file}的合并序列已保存至{new_combined_species_seqs}')
            print('\n')
        else:
            print(f'物种{file}没有可扩增的序列')

    # 处理CSV文件
    fp = pd.read_csv(input_path + '/sequence_info.csv', encoding='utf_8_sig')

    # 修改点3：只对遇到的GI设置0或1，未遇到的设为NaN
    fp['amplify'] = fp['GI'].apply(lambda x:
                                   0 if str(x) in file_seqs_needmove else
                                   1 if str(x) in encountered_gis else
                                   np.nan)

    seq_info_df = pd.DataFrame(seq_info)
    seq_info_df['GI'] = seq_info_df['GI'].astype('int64')

    merged_df = pd.merge(fp, seq_info_df, on='GI', how='left')
    merged_df.to_csv(output_path + '/ISPCR_info.csv', encoding='utf_8_sig', index=False)

    return print('finished ISPCR')


# ============================================================
# 测试优化后的函数
#
# 所有 primer 均按 5' -> 3' 输入。
# R 不要提前反向互补，因为程序内部已经转换为 reverse complement。
# ============================================================

# 12S-V5，修正后的方向
# isPCR_optimized(
#     input_path=r'D:\1\researches\By WJ\2中期\2024gz12stele02\big\Download',
#     F='TTAGATACCCCACTATGC',
#     R='TAGAACAGGCTCCTAG',
#     mis=1
# )

# Tele02 示例：
# isPCR_optimized(
#     input_path=r'F:\gdut\学习\3答辩\嵌入\dq\Download',
#     F='AAACTCGTGCCAGCCACC',
#     R='GGGTATCTAATCCCAGTTTG',
#     mis=1
# )

# Riaz-12S 示例：
# isPCR_optimized(
#     input_path=r'D:\1\researches\By WJ\2中期\2024gz12stele02\big\Download',
#     F='ACTGGGATTAGATACCCC',
#     R='TAGAACAGGCTCCTAG',
#     mis=1
# )
isPCR_optimized(
    input_path=r'F:\gdut\学习\3答辩\嵌入\dq\Download',
    F='GGWACWGGWTGAACWGTWTAYCCYCC',
    R='TAAACTTCAGGGTGACCAAAAAATCA',
    mis=1
)
