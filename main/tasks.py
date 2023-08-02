from datetime import datetime, timezone
from Bio import Entrez
from Bio import Medline
from celery import shared_task
import requests
import time
import json
import os
import math
import io

import oneai
import requests

from .serializers import *
from dm.settings import PARSER_EMAIL, MEDIA_URL, BERTTOPICAPI_URL, SUMMARISEAPI_URL
from Funcs import *


def get_path_to_file(pk, file_name):
    path_to_file = os.path.join('datasets', str(pk), file_name)
    print(path_to_file)
    if not os.path.exists(path_to_file):
        raise Exception('Not found current_file')

    return path_to_file

def check_permission(user_permission, IdList):
    allow_records = user_permission.all_records - user_permission.used_records
    print(f'Allow analise records = {allow_records}, Len recieved data = {len(IdList)}')
    if user_permission.used_records is not None:
        if not (len(IdList) > allow_records):
            user_permission.used_records += len(IdList)
        else:
            IdList = IdList[:allow_records]
            user_permission.used_records = user_permission.all_records
        user_permission.save()
    return IdList


@shared_task(bind=True)
def parse_records(self, query: str, count: int, new_task_id: int, permission_id: int, retmax: int = 10000):
    # Парсим полученные записи
    print("Start parsing...")
    print(f"Count={count}, Task id={new_task_id}")
    new_task = TaskSearch.objects.get(id=new_task_id)
    user_permission = UserPermissions.objects.get(id=permission_id)

    handle = Entrez.esearch(db="pubmed", sort='relevance', term=query, retmax=retmax)
    f = Entrez.read(handle)
    IdList = f['IdList']
    IdList = check_permission(user_permission, IdList)

    handle.close()

    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")

    records = []
    i = 0

    for record in Medline.parse(handle):
        data = parse_record(record)
        if data:
            records.append(data)

        if i % 1000 == 0:
            new_task.message = f'On {i}/{count} step...'
            new_task.save()

        i += 1

    new_task.save()
    print(len(records))
    data = {
        'search_ncbi': ArticleSerializer(records, many=True).data
    }
    handle.close()
    return data


# Удаление пробелов и Up-первого символа слов строки
def move_space(inp_str):

  new_str = inp_str.strip() # удаляем начальные и конечные пробелы из сроки
  while True:
     i = new_str.find(' ') # находим пробел
     if (i < 0):
      break
     a = new_str[:i]
     b = new_str[i+1:]
     b = b.capitalize()
     new_str = a + b

  return new_str


@shared_task(bind=True)
def analise_records(self, pk, params, new_task_id, user_permission_id):
    new_task = TaskAnalise.objects.get(id=new_task_id)
    new_task.message = 'Запущен процесс анализа, пожайлуста подождите...'
    new_task.save()

    response = requests.post(f'{BERTTOPICAPI_URL}/analise', json=params)
    print(response.status_code)
    if response.status_code != 200:
        raise Exception('Api is failed!')

    data = response.json()
    with open(get_path_to_file(pk, 'heapmap.json'), 'w') as f:
        json.dump(data['heapmap'], f)
    with open(get_path_to_file(pk, 'tematic_analise.json'), 'w') as f:
        user_permission = UserPermissions.objects.get(id=user_permission_id)
        check_permission(user_permission, data['tematic_analise'])
        json.dump(data['tematic_analise'], f)
    with open(get_path_to_file(pk, 'clust_graph.json'), 'w') as f:
        json.dump(data['clust_graph'], f)
    with open(get_path_to_file(pk, 'heirarchy.json'), 'w') as f:
        json.dump(data['heirarchy'], f)
    with open(get_path_to_file(pk, 'DTM.json'), 'w') as f:
        json.dump(data['DTM'], f)
    with open(get_path_to_file(pk, 'topics.json'), 'w') as f:
        json.dump(data['topics'], f)

    out_v = io.open(get_path_to_file(pk, 'vectors_OR.tsv'), 'w', encoding='utf-8')
    out_m = io.open(get_path_to_file(pk, 'metadata_OR.tsv'), 'w', encoding='utf-8')
    l = 0
    print('Create embeddings vectors tsv...')
    for emb in data['embeddings']:
        str_v = ''
        for e in emb:
            str_v += '\t' + str(e)
        str_v += "\n"
        out_v.write(str_v)
        l += 1

    l = 0
    print('Create embeddings metadata tsv...')
    for record in data['tematic_analise']:
        if (record['topic'] >= 0):
            tpcs = str(record['topic'])
        else:
            tpcs = '#'
        out_m.write(move_space(record['titl'][:90]) + '#T' + tpcs + '_D' + str(l + 1) + "\n")
        l += 1

    out_v.close()
    out_m.close()
    print('done!')

    return


@shared_task(bind=True)
def summarise_text(self, records):
    response = requests.post(f'{SUMMARISEAPI_URL}/summarise', json={'IdList': records})
    print(response.status_code)
    if response.status_code != 200:
        raise Exception('Api is failed!')

    sentences = ' '.join(response.json()['summary'])

    return sentences


@shared_task(bind=True)
def summarise_emb(self, pk, text):
    text = ' '.join([s for s in text])
    response = requests.post(f'{SUMMARISEAPI_URL}/summarise_emb', json={'text': text})
    print(response.status_code)
    if response.status_code != 200:
        raise Exception('Api is failed!')

    sentences = ' '.join(response.json()['summary'])

    return sentences


def current_record(record, **filters):
    if 'score' in filters:
        if record['score'] < filters['score']:
            return False
    return True


def filter_record(record, **filters):
    if 'date' in filters:
        if filters['date']:
            period = datetime.now().year - int(record.pdat.split(' ')[0])
            print(f"{period} ?> {filters['date']}")
            if period > filters['date']:
                return False
    if not ('type' in filters and len(filters['type']) > 0):
        return True
    else:
        res = False
        current_types = record.pt.split('; ')
        print(current_types, filters['type'])
        for filter_type in filters['type']:
            if filter_type in current_types:
                res = True
    return res


@shared_task(bind=True)
def get_ddi_articles(self, query, new_task_id, user_permission_id, **kwargs):

    url = f'https://www.ncbi.nlm.nih.gov/research/litsense-api/api/?query={query}&rerank=true'
    r = requests.get(url)
    print(r.status_code)
    if (r.status_code > 299 or r.status_code < 200):
        raise ConnectionError(f'Подключение к эмбедингам ncbi прошло неудачно с {r.status_code} ошибкой')

    user_permission = UserPermissions.objects.get(id=user_permission_id)
    records = check_permission(user_permission, json.loads(r.text))
    records_id = {record['pmid']: [round(record['score'], 2), record['section'], record['text']] for record in records if current_record(record, **kwargs)}

    handle = Entrez.efetch(db="pubmed", id=[k for k in records_id], rettype="medline", retmode="text")

    records = []
    len_records = len(records_id)
    new_task = TaskAnalise.objects.get(id=new_task_id)
    i = 0

    for record in Medline.parse(handle):
        data = parse_record(record)
        if not filter_record(data, **kwargs):
            i += 1
            continue
        if data:
            record = ArticleSerializer(data, many=False).data
            record['annotations'] = None
            record['query_number'] = kwargs['number_of_query']
            record['score'], record['section'], record['text'] = records_id[record['uid']]

            records.append(record)

        if i % 10 == 0:
            new_task.message = f'On {i}/{len_records} step...'
            new_task.save()

        i += 1

    print(len(records))
    handle.close()
    return records

@shared_task(bind=True)
def markup_artcile(self, record):
    annotations = get_pmid(record['uid'])
    if annotations is None:
        annotations = get_annotations(record['tiab'])

    if annotations is None:
        record['annotations'] = None
    else:
        for annotation in annotations['annotations']:
            if math.isnan(annotation['prob']):
                annotation['prob'] = None
    record['tiab'] = annotations['text']
    record['annotations'] = annotations['annotations']

    return record

@shared_task(bind=True)
def plot_graph_associations(self, IdList, pk, max_size=200):
    # Start plot graph on tematic analise

    if len(IdList) > max_size:
        IdList = IdList[:max_size - 1]

    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")
    articles = []

    for record in Medline.parse(handle):
        data = parse_record(record)
        if data:
            articles.append(data)

    handle.close()

    data = get_uniq_info_for_authors(articles)
    f = open(get_path_to_file(pk, 'info_graph.json'), 'w')
    json.dump(data, f)
    f.close()
    data = get_uniq_info_for_institutes(articles)
    f = open(get_path_to_file(pk, 'info_graph_affiliations.json'), 'w')
    json.dump(data, f)
    f.close()
    data = get_uniq_info_for_countries(articles)
    f = open(get_path_to_file(pk, 'info_graph_countries.json'), 'w')
    json.dump(data, f)
    f.close()
    data = get_uniq_info_for_other(articles)
    f = open(get_path_to_file(pk, 'info_graph_journals.json'), 'w')
    json.dump(data, f)
    f.close()
    return


