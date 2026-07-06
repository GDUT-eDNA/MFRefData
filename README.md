# 📦 Project Name

> MFRefData

---

本说明文档提供英文和中文两个版本

点击链接前往不同版本[English](#english-version) | [中文](#中文版本)  

---  

<a id="english-version"></a>  

## English Version

## 📖 1. Project Introduction

This project is designed for:

- **Purpose:** A Supervised and Visualization-Enabled Toolkit for Reducing Reference Database Uncertainty in eDNA Biomonitoring.

- **Features:**
  
  - Standardized and workflow-based
  
  - Visualized and supervised
  
  - Compatible with multiple formats and extensible

---

## ⚙️ 2. Environment and Dependencies

#### Runtime Environment:

- Python 3.x

- OS: Windows

#### Dependencies:

```python
pandas,beautifulsoup4,lxml,biopython,
matplotlib,squarify,matplotlib-venn,requests, etc.
```

### Environment Setup:

```python
# Clone the project repository
git clone https://github.com/GDUT-eDNA/MFRefData.git
# Create a Conda environment
conda create -n MFRefData python=3.10 -y
# Install dependencies
pip install pandas beautifulsoup4 lxml requests
pip install biopython matplotlib squarify matplotlib-venn
```

---

## ⚙️ 3. Function Description

```text
project/
    │── taxdmp/             # Folder containing files required for annotation
    │── vsearch/            # Folder containing files required for clustering
    │── cluster.py          # Clustering
    │── db.py               # Sequence downloading
    │── extend.py           # Extension
    │── filter.py           # Filtering
    │── output.py           # Output
    │── remove_primers.py   # Primer removal
    │── totax.py            # Annotation
    │── visual_bar.py       # Visualization of species composition and distribution
    │── visual_pie.py       # Visualization of amplification results
    │── visual_venn.py      # Visualization of database comparisons
```

---

## ⚙️ 4. Quick Start

#### 4.1 Data Preparation

- Species list
- Create a virtual environment

#### 4.2 Example Usage

Run the following command:

```bash
python main.py -list example.txt -region 12s -type qiime2
```

#### 4.3 Quick Usage

1. Place the `species list.txt` file in the `work` folder.

2. Run the following command:

```bash
python main.py -list <species_list_file> -region <amplification_region> -type <output_type>
```

**Parameter Description:**

`-list`  
The species list used to download the required database. Database sequences will be downloaded according to this species list. The file should be placed in the `work` folder.

`-region`  
The target region for database downloading. Available options include `12s`, `16s`, `18s`, `coi`, etc.

`-type`  
The output database format. Available options include `qiime2`, `sintax`, `RDP`, `BLAST`, `kraken2`, `DADA2`, and `mothur`.

**Examples:**

```bash
python main.py -list species1.txt -region 12s -type sintax
python main.py -list species2.txt -region 18s -type qiime2
```

### 4.4 Output Files (Using QIIME 2 as an Example)

1. **Reference Sequence Database (FASTA format):** `output_qiime2.fasta`
   
   1. This file contains all reference nucleotide sequences used for taxonomic classification and sequence alignment.
   
   2. The header of each sequence begins with a unique identifier, such as `>ASV_1` or `>Gene_XXX`, which serves as the primary key for that sequence.

2. **Taxonomic Lineage Mapping File (TXT format):** `output_qiime2.txt`
   
   1. This file is a tab-delimited plain-text table.
   
   2. It contains two columns:
      
      1. The first column contains the sequence identifier, which must correspond exactly to the sequence ID in the FASTA file.
      
      2. The second column contains the complete taxonomic lineage of the sequence. Taxonomic ranks are arranged in the following order: kingdom, phylum, class, order, family, genus, and species. Each rank is separated by a semicolon (`;`), for example:
      
      `k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria`
      
      Missing taxonomic ranks should be represented by `__` or `unclassified`.

---

<a id="中文版本"></a>

## 中文版本

## 📖 1. 项目介绍

此项目用于：

- 目的: 一种用于降低环境DNA生物监测中参考数据库不确定性的监督式可视化工具包
- 特点:
  - 标准化、流程化
  - 可视化、可监督 
  - 多格式兼容、可扩展

---

## ⚙️ 2. 环境与依赖

#### 运行环境：

- Python 3.x
- OS: Windows 

#### 依赖：

```python
pandas、beautifulsoup4、lxml、biopython、
matplotlib、squarify、matplotlib-venn、requests等
```

#### 环境建立：

```python
#项目下载
git clone https://github.com/GDUT-eDNA/MFRefData.git
#创建环境
conda create -n MFRefData python=3.10 -y
#安装依赖
pip install pandas beautifulsoup4 lxml requests
pip install biopython matplotlib squarify matplotlib-venn
```

----

## ⚙️3.函数介绍

```text
project/  
        │── taxdmp/           #注释所需文件夹
        │── vsearch/           #聚类所需文件夹
        │── cluster.py         #聚类 
        │── db.py                #下载序列 
        │── extend.py        #扩展 
        │── filter.py            #过滤
        │── output.py        #输出
        │── remove_primers.py #去引物 
        │── totax.py           #注释
        │── visual_bar.py          #物种组成与分布可视化
        │── visual_pie.py          #扩增情况可视化
        │── visual_venn.py          #数据库对比可视化
```

---

## ⚙️4.快速使用

#### 4.1 数据准备：

- 物种名录
- 创建虚拟环境

#### 4.2 示例使用：

输入指令：

```bash
python main.py -list example.txt -region 12s -type qiime2
```

#### 4.3 快速使用：

1. 将`物种名录.txt`放入`work`文件夹

2. 输入指令：

```bash
python main.py -list <物种名录文件> -region <扩增区域> -type <输出类型>
```

**参数解释**：

`-list`    为所需下载数据库的物种名录，根据此物种名录进行数据库下载（放入work文件夹）

`-region`    为目标数据库下载区域，可选有`12s` `16s` `18s` `coi`等

`-type`    为输出数据库格式，可选有`qiime2` `sintax` `RDP` `BLAST` `kraken2` `DADA2` `mothur`

**示例:**

```bash
python main.py -list species1.txt -region 12s -type sintax
python main.py -list species2.txt -region 18s -type qiime2
```

#### 4.4 输出：（以qiime2为例）

1. **参考序列数据库（FASTA 格式）** ：即 `output_qiime2.fasta` 文件
   
   1. 存储用于分类比对的所有参考核苷酸序列
   
   2. 每条序列的头部（Header）均以唯一标识符（如 `>ASV_1` 或 `>Gene_XXX`）作为序列主键

2. **分类谱系映射文件（TXT 文本格式）**：即 `output_qiime2.txt` 文件
   
   1. 为制表符分隔的纯文本表格
   
   2. 该文件包含两列：
      
      1. 第一列为序列标识符（与 FASTA 文件中的 ID 严格一一对应）
      
      2. 第二列为该序列完整的分类学谱系。谱系按 界、门、纲、目、科、属、种的层级顺序排列，各层级间以分号（`;`）分隔（例如：`k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria`），缺失的层级以 `__` 或 `unclassified` 占位
