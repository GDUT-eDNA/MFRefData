import requests
import pandas as pd
import os
import time
from bs4 import BeautifulSoup
import re
from Bio import SeqIO
from io import StringIO

tool = 'MFdb'
geo_pattern = re.compile(r'/geo_loc_name\s*=\s*"([^"]+)"')
reference_pattern = re.compile('REFERENCE')
id_pattern = re.compile(r'VERSION\s+(\S+)')


def download_NCBI(path, region, count):
    filepath = path
    NCBI_search = region
    download_count_per_species = count
    output_dir = os.path.dirname(path)
    output_file = os.path.join(output_dir, 'Download')
    output_dir_sequence = os.path.join(output_file, 'sequences')

    if not os.path.exists(path):
        print(f"{path} is not existed, please check")
        return

    os.makedirs(output_file, exist_ok=True)
    os.makedirs(output_dir_sequence, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        "User-Agent": "MFdb/2.0,2532707857@qq.com",
        "Accept": "text/html,application/xhtml+xml",
        "Accept-Encoding": "gzip, deflate",
        "Connection": "keep-alive"
    })

    # 读取物种列表
    print('reloading species file')
    species_df = pd.read_csv(filepath, sep="\t", header=None)
    species_list = species_df[0].tolist()
    print(f"There are {len(species_list)} species in the list {path}:")
    print(f"{species_list} \n")

    taxid_cache = {}
    start = time.time()


    def get_taxid(species_name):
        if species_name in taxid_cache:
            return taxid_cache[species_name]

        params = {
            "db": "taxonomy",
            "term": f'"{species_name}"[Organism]',
            "retmode": "xml"
        }

        try:
            request_start = time.time()
            response = session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params=params,
                timeout=(10, 30)
            )
            response.raise_for_status()

            soup = BeautifulSoup(response.text, "xml")
            taxid = soup.find("Id").text if soup.find("Id") else None

            taxid_cache[species_name] = taxid

            # 保持NCBI速率限制：3次/秒
            request_end = time.time() - request_start
            print('taxid:',request_end)
            if request_end < 0.34:
                time.sleep(0.34 - request_end)
                final_time = time.time() - request_start
                print('taxid final',final_time)
            time.sleep(0.07)
            return taxid
        except Exception as e:
            print(f"get {species_name} taxid is fall: {str(e)}")
            return None

    sequence_info_all = []
    num_seq = 0

    for species in species_list:
        # 获取TaxID
        taxid = get_taxid(species)
        if taxid is None:
            print(f'{species} has no taxid, please check')
            print('\n')
            continue

        # 搜索序列
        search_term = f'"{species}"[Organism] AND {NCBI_search}'
        search_params = {
            "db": "nuccore",
            "term": search_term,
            "retmax": download_count_per_species,
            "retmode": "xml",
            "usehistory": "y",
        }

        try:
            response = session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params=search_params,
                stream=True,
                timeout=(10, 30)
            )
            request_start = time.time()
            response.raise_for_status()
            soup = BeautifulSoup(response.text, "xml")

            count_element = soup.find('Count')
            webenv_element = soup.find('WebEnv').text
            query_element = soup.find('QueryKey').text

            seq_exist = int(count_element.text)

            if seq_exist < download_count_per_species and seq_exist != 0:
                print(
                    f"{species} in {search_term} has {seq_exist} sequences, can not download enough {download_count_per_species}.\n{seq_exist} sequences will be download")

            elif seq_exist > download_count_per_species:
                print(
                    f"{species} in {search_term} has enough {seq_exist} sequences, {download_count_per_species} download enough")

            elif seq_exist == 0:
                print('\n')
                request_end = time.time() - request_start
                if request_end < 0.34:
                    time.sleep(0.34 - request_end)
                continue

            real_download_num = min(seq_exist,download_count_per_species)
            taxid_path = os.path.join(output_dir_sequence, f"{species}_{taxid}")
            os.makedirs(taxid_path, exist_ok=True)


            GI_accessions = [id_tag.text for id_tag in soup.find_all("Id")]
            num_seq = num_seq + len(GI_accessions)

            # 保持NCBI速率限制
            request_end = time.time() - request_start
            if request_end < 0.34:
               time.sleep(0.34 - request_end)

            # 下载序列
            download_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                "db": "nuccore",
                "rettype": "gbwithparts",
                "retmode": "text",
                "query_key": query_element,
                "WebEnv": webenv_element,
                "retstart":0,
                "retmax":real_download_num
            }

            content = ''
            request_start = time.time()

            try:
                response = session.get(download_url, params=params,stream=True, timeout=(12, 20))
                response.raise_for_status()
                content_chunks = []
                for chunk in response.iter_content(chunk_size=131072, decode_unicode=True):
                    content_chunks.append(chunk)

                content = ''.join(content_chunks)

            except Exception as e:
                request_end = time.time()
                first_attempt = request_end - request_start
                print(e)
                print(f'try download {taxid} again with {len(GI_accessions)}:')
                print(f'{GI_accessions}')
                print("first:",first_attempt)

                part_real_download_num = real_download_num // 2
                part_start = 0

                content_part = []
                for j in range(2):
                    try:
                        params = {
                            "db": "nuccore",
                            "rettype": "gbwithparts",
                            "retmode": "text",
                            "query_key": query_element,
                            "WebEnv": webenv_element,
                            "retstart": part_start,
                            "retmax": part_real_download_num
                        }
                        part_start = part_real_download_num + part_start
                        part_real_download_num = real_download_num - part_real_download_num

                        request_start = time.time()
                        response = session.get(download_url, params=params,stream=True, timeout=(12, 20))
                        response.raise_for_status()

                        content_chunks = []
                        for chunk in response.iter_content(chunk_size=131072, decode_unicode=True):
                            content_chunks.append(chunk)

                        content_part.append(''.join(content_chunks))

                        request_end = time.time() - request_start
                        if request_end < 0.34:
                            time.sleep(0.34 - request_end)

                    except Exception as e:
                        print(e)
                        break

                if len(content_part) == 2:
                    content = ''.join(content_part)
                    print('second success!')

            if content == '':
                print(f'{taxid} download finally failed')
                print('\n')
                continue
            else:
                records = content.rstrip().split('\n//\n')
                # 处理每个记录
                if len(records) != len(GI_accessions):
                    print('Index out of range')
                    continue
                else:
                    for i, record in enumerate(records):
                        # 使用预编译的正则表达式
                        geo_match = geo_pattern.search(record)
                        geo = geo_match.group(1).replace('\n', ' ').replace('\r', ' ').strip() if geo_match else 'na'
                        geo = re.sub(r'\s+', '', geo)
                        ref_count = len(reference_pattern.findall(record))
                        id_match = id_pattern.search(record)
                        sequence_id = id_match.group(1) if id_match else 'na'

                        GI = GI_accessions[i]
                        sequence_info_all.append({
                            'scientific name': species,
                            'taxid': taxid,
                            'name': sequence_id,
                            'geography': geo,
                            'reference': ref_count,
                            'GI': GI
                        })

                        try:
                            record = record.rstrip() + '\n//\n'
                            for seq_record in SeqIO.parse(StringIO(record),format='gb'):
                                fasta_filename = os.path.join(taxid_path, f'{GI}.fasta')
                                with open(fasta_filename, 'w', encoding='utf-8') as f:
                                    SeqIO.write(seq_record,f,format='fasta')

                        except Exception as e:
                            print(f"Failed to parse sequence {GI}: {str(e)}")
                            continue

                # 保持NCBI速率限制
                request_end = time.time() - request_start
                print('download:',request_end)
                print('\n')
                if request_end < 0.34:
                    time.sleep(0.34 - request_end)
        except Exception as e:
            print(f" {species} need be focused: {str(e)}")
            print('\n')
            continue


    # 保存CSV
    output_path = os.path.join(output_file, "sequence_info.csv")
    df = pd.DataFrame(sequence_info_all)
    df.to_csv(output_path, index=False, encoding="utf-8-sig")

    end = time.time() - start
    print(f"\nTotal time: {end:.2f} seconds")
    print(f"Total sequences: {num_seq}")

    return 0

if __name__ == "__main__":
    work_dir = './work'
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    Species_list = example.txt
    Species_region = '12s'
    download_NCBI(
        path=os.path.join(work_dir, Species_list),
        region=Species_region,
        count=100
    )
