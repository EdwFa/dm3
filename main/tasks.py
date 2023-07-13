from Bio import Entrez
from Bio import Medline
from celery import shared_task
import time
import json
import os
import math

import oneai
import requests

from django_celery_results.models import TaskResult

from .serializers import *
from dm.settings import PARSER_EMAIL, MEDIA_URL
from Funcs import *


def check_working_task(request, Task, **kwargs):
    tasks = Task.objects.filter(**kwargs)
    for task in tasks:
        try:
            worker = TaskResult.objects.get(task_id=task.task_id)
        except:
            task.delete()
            continue
        else:
            if worker.status == 'PROGRESS' or worker.status == 'STARTED':
                return True

    return False


def get_path_to_file(username, file_name):
    path_to_file = os.path.join('datasets', username, file_name)
    print(path_to_file)
    if not os.path.exists(path_to_file):
        raise Exception('Not found current_file')

    return path_to_file


@shared_task(bind=True)
def parse_records(self, query: str, count: int, new_task_id: int, retmax: int = 10000):
    # Парсим полученные записи
    print("Start parsing...")
    print(f"Count={count}, Task id={new_task_id}")
    new_task = TaskSearch.objects.get(id=new_task_id)
    Entrez.email = PARSER_EMAIL  # Говорю NCBI кто я есть

    handle = Entrez.esearch(db="pubmed", sort='relevance', term=query, retmax=retmax)
    f = Entrez.read(handle)
    IdList = f['IdList']
    handle.close()

    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")

    records = []
    i = 0

    for record in Medline.parse(handle):
        data = parse_record(record)
        if data:
            records.append(data)

        if i % 100 == 0:
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

def create_analise_task(request, **kwargs):
    if check_working_task(request, TaskAnalise, user=kwargs['user'], status=0):
        return None

    task = TaskAnalise.objects.create(**kwargs)
    return task


def check_working_analise_task(request, type_status):
    if check_working_task(request, TaskAnalise, status=0, user=request.user, type_analise=type_status):
        return True
    # if tasks.count() != 0:
    #     return True

    return False

@shared_task(bind=True)
def analise_records(self, username, IdList, new_task_id):
    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")
    records = ArticleSerializer([parse_record(record) for record in Medline.parse(handle) if record], many=True).data
    handle.close()


    new_task = TaskAnalise.objects.get(id=new_task_id)
    new_task.message = 'Предобработка данных перед передачей модели...'
    new_task.save()

    articles = create_clear_articles(records)

    new_task.message = 'Анализ записей моделью...'
    new_task.save()
    topics, props, embeddings = analise_articles(articles)
    records = return_results(records, topics, props)

    new_task.message = 'Вывод результатов и отрисовка графиков...'
    new_task.save()
    graph = return_clust_graph([rec['titl'] for rec in records], embeddings)
    count_topics = len(set(topics))
    print("Count topic = ", count_topics)
    n_clusters = 10
    if count_topics > 1:
        if n_clusters >= count_topics - 2:
            n_clusters = int(count_topics / 2)

        print(n_clusters)
        heapmap = return_heapmap(n_clusters=n_clusters)
        with open(get_path_to_file(username, 'heapmap.json'), 'w') as f:
            f.write(heapmap)
    else:
        with open(get_path_to_file(username, 'heapmap.json'), 'w') as f:
            json.dump(None, f)

    with open(get_path_to_file(username, 'tematic_analise.json'), 'w') as f:
        json.dump(records, f)
    with open(get_path_to_file(username, 'clust_graph.json'), 'w') as f:
        f.write(graph)
    with open(get_path_to_file(username, 'heirarchy.json'), 'w') as f:
        try:
            f.write(return_heirarchy())
        except:
            json.dump(None, f)
    return None

@shared_task(bind=True)
def summarise_text(self, records):
    text = ' '.join([rec['titl'] + rec['tiab'] for rec in records])
    print(len(text), text)
    api_key = "0f4f9d19-644d-4577-94af-48abb689be60"
    oneai.api_key = api_key
    pipeline = oneai.Pipeline(steps=[
        oneai.skills.Summarize(min_length=20),
    ])

    output = pipeline.run(text)
    return output.summary.text

@shared_task(bind=True)
def summarise_emb(self, username):
    with open(get_path_to_file(username, 'embeddings.json'), 'r') as f:
        data = json.load(f)['data']
    text = ' '.join([rec['text'] for rec in data])
    print(len(text), text)
    api_key = "0f4f9d19-644d-4577-94af-48abb689be60"
    oneai.api_key = api_key
    pipeline = oneai.Pipeline(steps=[
        oneai.skills.Summarize(min_length=20),
    ])

    output = pipeline.run(text)
    return output.summary.text


def current_record(record, **filters):
    if 'score' in filters:
        if record['score'] < filters['score']:
            return False
    return True


@shared_task(bind=True)
def get_ddi_articles(self, query, new_task_id, **kwargs):

    url = f'https://www.ncbi.nlm.nih.gov/research/litsense-api/api/?query={query}&rerank=true'
    r = requests.get(url)
    records_id = {record['pmid']: [round(record['score'], 2), record['section'], record['text']] for record in json.loads(r.text) if current_record(record, **kwargs)}

    handle = Entrez.efetch(db="pubmed", id=[k for k in records_id], rettype="medline", retmode="text")

    records = []
    len_records = len(records_id)
    new_task = TaskAnalise.objects.get(id=new_task_id)
    i = 0

    for record in Medline.parse(handle):
        data = parse_record(record)
        if data:
            record = ArticleSerializer(data, many=False).data
            annotations = get_pmid(record['uid'])
            if annotations is None:
                annotations = get_annotations(record['tiab'])

            print(annotations is None)
            if annotations is None:
                record['annotations'] = None
                continue
            else:
                for annotation in annotations['annotations']:
                    if math.isnan(annotation['prob']):
                        annotation['prob'] = None
            record['tiab'] = annotations['text']
            record['annotations'] = annotations['annotations']
            record['score'], record['section'], record['text'] = records_id[record['uid']]

            records.append(record)

        if i % 10 == 0:
            new_task.message = f'On {i}/{len_records} step...'
            new_task.save()

        i += 1

    print(len(records))
    data = {
        'embeddings': records
    }
    handle.close()
    return data



