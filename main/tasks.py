from Bio import Entrez
from Bio import Medline
from celery import shared_task
import time
import json
import oneai
import requests


from .serializers import *
from dm.settings import PARSER_EMAIL, MEDIA_URL
from Funcs import *


def create_task(**kwargs):
    task = Task.objects.create(**kwargs)
    return task


def check_working_task(request):
    tasks = Task.objects.filter(status=0, user=request.user)
    if tasks.count() != 0:
        return True

    return False


@shared_task(bind=True)
def parse_records(self, query: str, count: int, new_task_id: int, retmax: int = 10000):
    # Парсим полученные записи
    print("Start parsing...")
    print(f"Count={count}, Task id={new_task_id}")
    new_task = Task.objects.get(id=new_task_id)
    Entrez.email = PARSER_EMAIL  # Говорю NCBI кто я есть

    handle = Entrez.esearch(db="pubmed", sort='relevance', term=query, retmax=retmax)
    f = Entrez.read(handle)
    IdList = f['IdList']
    handle.close()

    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")

    records = []

    for record in Medline.parse(handle):
        data = parse_record(record)
        if data:
            records.append(data)
            new_task.save()

    print(len(records))

    handle.close()
    data = {
        'data': ArticleSerializer(records, many=True).data,
        'message': 'no one task in progress'
    }
    return data

def create_analise_task(**kwargs):
    if TaskAnalise.objects.filter(user=kwargs['user'], status=0).count() != 0:
        return None
    task = TaskAnalise.objects.create(**kwargs)
    return task


def check_working_analise_task(request):
    tasks = TaskAnalise.objects.filter(status=0, user=request.user, type_analise=0)
    if tasks.count() != 0:
        return True

    return False

@shared_task(bind=True)
def analise_records(self, records):
    articles = create_clear_articles(records)
    topics, props, embeddings = analise_articles(articles)
    records = return_results(records, topics, props)
    graph = return_clust_graph([rec['titl'] for rec in records], embeddings)
    count_topics = len(set(topics))
    n_clusters = 10
    if n_clusters > count_topics:
        n_clusters = count_topics / 2
    heapmap = return_heapmap(n_clusters=n_clusters)
    heirarchy = return_heirarchy()

    with open('test_analise_json.json', 'w') as f:
        json.dump(records, f)
    with open('test_clust_graph.json', 'w') as f:
        f.write(graph)
    with open('test_heapmap.json', 'w') as f:
        f.write(heapmap)
    with open('test_heirarchy.json', 'w') as f:
        f.write(heirarchy)
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
def get_ddi_articles(self, query):
    url = f'https://www.ncbi.nlm.nih.gov/research/litsense-api/api/?query={query}&rerank=true'
    r = requests.get(url)
    jsn = json.loads(r.text)
    return jsn



