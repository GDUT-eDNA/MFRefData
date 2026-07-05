# 📦 Project Name

> MFRefData

---
This README is available in English and Chinese. Click the link to access different versions.

本说明文档提供英文和中文两个版本，点击链接前往不同版本[English](#english-version) | [中文](#中文版本)  

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
pandas, beautifulsoup4, lxml, etc.
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
pandas、beautifulsoup4、lxml等
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
