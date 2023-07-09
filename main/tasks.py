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
    i = 0

    for record in Medline.parse(handle):
        data = parse_record(record)
        if data:
            records.append(data)

        if i % 10 == 0:
            new_task.message = f'On {i}/{count} step...'
            new_task.save()

        i += 1

    new_task.save()

    print(len(records))

    handle.close()
    data = {
        'data': ArticleSerializer(records, many=True).data,
        'message': 'Done!'
    }
    return data

def create_analise_task(**kwargs):
    if TaskAnalise.objects.filter(user=kwargs['user'], status=0).count() != 0:
        return None
    task = TaskAnalise.objects.create(**kwargs)
    return task


def check_working_analise_task(request, type_status):
    tasks = TaskAnalise.objects.filter(status=0, user=request.user, type_analise=type_status)
    if tasks.count() != 0:
        return True

    return False

@shared_task(bind=True)
def analise_records(self, IdList, new_task_id):
    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")
    records = ArticleSerializer([parse_record(record) for record in Medline.parse(handle) if record], many=True).data
    handle.close()


    new_task = TaskAnalise.objects.get(id=new_task_id)
    new_task.message = 'Start clearing articles for model...'
    new_task.save()

    articles = create_clear_articles(records)

    new_task.message = 'Analise articles in model...'
    new_task.save()
    topics, props, embeddings = analise_articles(articles)
    records = return_results(records, topics, props)

    new_task.message = 'Printing graphics after analising...'
    new_task.save()
    graph = return_clust_graph([rec['titl'] for rec in records], embeddings)
    count_topics = len(set(topics))
    print("COunt topic = ", count_topics)
    n_clusters = 10
    if count_topics > 1:
        if n_clusters > count_topics:
            n_clusters = count_topics - 1
        heapmap = return_heapmap(n_clusters=n_clusters)
        with open('test_heapmap.json', 'w') as f:
            f.write(heapmap)
    else:
        with open('test_heapmap.json', 'w') as f:
            json.dump(None, f)

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
def get_ddi_articles(self, query, new_task_id):
    # # оnly for test check status
    # new_task = TaskAnalise.objects.get(id=new_task_id)
    # for i in range(10):
    #     new_task.message = f'On {i}/10 step...'
    #     new_task.save()
    #     time.sleep(1)
    # return None

    url = f'https://www.ncbi.nlm.nih.gov/research/litsense-api/api/?query={query}&rerank=true'
    r = requests.get(url)
    records_id = {record['pmid']: [round(record['score'], 2), record['section'], record['text']] for record in json.loads(r.text)}

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
                if annotations is None:
                    continue
            record['tiab'] = markup_text(**annotations)
            record['score'], record['section'], record['text'] = records_id[record['uid']]

            records.append(record)

        if i % 10 == 0:
            new_task.message = f'On {i}/{len_records} step...'
            new_task.save()

        i += 1

    print(len(records))

    handle.close()
    return records



