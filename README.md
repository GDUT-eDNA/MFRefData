# 📦 项目名称

> MFRefData

---

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
