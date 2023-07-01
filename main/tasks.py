from django.core.exceptions import ObjectDoesNotExist

from Bio import Entrez
from Bio import Medline
from datetime import date
from celery import shared_task
import numpy
import time

from .models import Article, Task
from dm.settings import PARSER_EMAIL, MEDIA_URL


# Запрос для начального запуска как пример
default_query = "{0} and {1}:{2}[dp] and Clinical Trial[pt]".format('covid-19', '2000/01/01',
                                                                    date.today().strftime("%Y/%m/%d"))

def create_query(request):
    # Преобразуем наши ключи в запрос pubmed
    print([i for i in request.GET.items()])

    if len([i for i in request.GET.items()]) == 0:
        return default_query

    # Содаем запрос с текстом запроса
    query = ""
    text = request.GET.get('search_field', None)
    if text:
        query = f'{text}'

    # Добавляем даты
    start_date = request.GET.get('dateStart', "1900-01-01")
    end_date = request.GET.get('dateStop', date.today().strftime("%Y-%m-%d"))
    query = f"{query} and {start_date.replace('-', '/')}:{end_date.replace('-', '/')}[dp]"

    # Добавляем гендеры
    genders = request.GET.getlist('Gender', None)
    if genders:
        query_genders = " OR ".join([f'{i}[mh]' for i in genders])
        query = f"{query} AND ({query_genders})"

    # Добавляем типы документов
    types = request.GET.getlist('Type', None)
    if types:
        query_types = " OR ".join([f'{i}[pt]' for i in types])
        query = f"{query} AND ({query_types})"

    # Добавляем Возрастные группы
    olds = request.GET.getlist('Old', None)
    if olds:
        query_olds = " OR ".join(olds)
        query = f"{query} AND ({query_olds})"

    return query

def get_records(query: str):
    # Получаем наши записи
    Entrez.email = "e.p@d_health.pro"  # Говорю NCBI кто я есть
    handle = Entrez.esearch(db="pubmed", sort='relevance', term=query)
    record = Entrez.read(handle)

    handle.close()

    str_query = record['QueryTranslation']

    try:
        translation_stack = record['TranslationStack']
    except:
        translation_stack = 'null/null/null'

    records_count = int(record['Count'])

    print(f"All records: {records_count}")

    return str_query, translation_stack, records_count

def create_task(**kwargs):
    task = Task.objects.create(**kwargs)
    return task

def check_working_task(request):
    tasks = Task.objects.filter(status=0, user=request.user)
    if tasks.count() != 0:
        return True

    return False

def parse_record(record):
    # Парсим полученный словарь записи
    if not ('PMID' in record): # Проверяем на наличие pmid если его нет не сохраняем данные
        return

    try:
        return Article.objects.get(uid=record['PMID'])
    except ObjectDoesNotExist:
        pass

    data = Article()
    data.uid = record['PMID']
    data.url = 'https://pubmed.ncbi.nlm.nih.gov/' + record['PMID'] + '/'
    if 'TI' in record:
        data.titl = record['TI']
    else:
        data.titl = "без title"

    if 'AID' in record:
        data.aid = ' '.join([''.join(i.split(' ')[:-1]) for i in record['AID']])
        data.urlaid = ' '.join([f'https://sci-hub.do/{aid}' for aid in data.aid.split(' ')])
    else:
        data.aid = "Без AID"
        data.urlaid = "Без url AID"

    if 'JT' in record:
        data.jour = record['JT']
    else:
        data.jour = "без JT"
    if 'DP' in record:
        data.pdat = record['DP']
    else:
        data.pdat = "без DP"
    if 'FAU' in record:
        data.auth = '; '.join(record['FAU'])
    else:
        data.auth = "не указаны"

    if 'AD' in record:
        data.affl = ';;; '.join(record['AD'])
    else:
        data.affl = "без AD"

    if 'PT' in record:
        data.pt = '; '.join(record['PT'])
    else:
        data.pt = "без PT"

    if 'PL' in record:
        data.pl = record['PL']
    else:
        data.pl = "без PL"

    if 'MH' in record:
        data.mesh = ';;; '.join(record['MH'])
    else:
        data.mesh = "без MH"

    if 'AB' in record:
        data.tiab = record['AB']
    else:
        data.tiab = "без AB"

    if 'OT' in record:
        data.ot = ';;; '.join(record['OT'])
    else:
        data.ot = "без OT"

    if 'RN' in record:
        data.rn = ';;; '.join(record['RN'])
    else:
        data.rn = "без RN"

    data.save()
    return data

# def parse_record(record):
#     # Парсим полученный словарь записи
#     if not ('PMID' in record): # Проверяем на наличие pmid если его нет не сохраняем данные
#         return
#
#     try:
#         return Article.objects.get(uid=record['PMID'])
#     except ObjectDoesNotExist:
#         pass
#
#     data = Article()
#     data.uid = record['PMID']
#     data.url = 'https://pubmed.ncbi.nlm.nih.gov/' + record['PMID'] + '/'
#     try:
#         data.titl = record['TI']
#     except:
#         data.titl = "без title"
#
#     try:
#         data.aid = ' '.join([''.join(i.split(' ')[:-1]) for i in record['AID']])
#         data.urlaid = ' '.join([f'https://sci-hub.do/{aid}' for aid in data.aid.split(' ')])
#     except:
#         data.aid = "Без AID"
#         data.urlaid = "Без url AID"
#
#     try:
#         data.jour = record['JT']
#     except:
#         data.jour = "без JT"
#     try:
#         data.pdat = record['DP']
#     except:
#         data.pdat = "без DP"
#     try:
#         data.auth = '; '.join(record['FAU'])
#     except:
#         data.auth = "не указаны"
#
#     try:
#         data.affl = ';;; '.join(record['AD'])
#     except:
#         data.affl = "без AD"
#
#     try:
#         data.pt = '; '.join(record['PT'])
#     except:
#         data.pt = "без PT"
#
#     try:
#         data.pl = record['PL']
#     except:
#         data.pl = "без PL"
#
#     try:
#         data.mesh = ';;; '.join(record['MH'])
#     except:
#         data.mesh = "без MH"
#
#     try:
#         data.tiab = record['AB']
#     except:
#         data.tiab = "без AB"
#
#     try:
#         data.ot = ';;; '.join(record['OT'])
#     except:
#         data.ot = "без OT"
#
#     try:
#         data.rn = ';;; '.join(record['RN'])
#     except:
#         data.rn = "без RN"
#
#     data.save()
#     return data

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

    print(len(records))
    new_task.articles.add(*records)

    handle.close()


def get_uniq_info_for_graph(articles, author_on_article=10):
    authors = []
    print('Get authors...')
    for article in articles:
        article_authors = article.auth.split('; ')
        if len(article_authors) > author_on_article:
            article_authors = article_authors[:author_on_article]
        for aa in article_authors:
            authors.append(aa)

    authors = set(authors)
    size = len(authors)
    print('Count of authors = ', size)
    print('Create nodes...')
    nodes = {k: i for i, k in enumerate(authors)}

    print('Getting info of edges...')
    edges = numpy.zeros((size, size), dtype='uint8')
    for article in articles:
        authors = article.auth.split('; ')
        if len(authors) > author_on_article:
            authors = authors[:author_on_article]
        for author in authors:
            for a in authors:
                edges[nodes[author]][nodes[a]] += 1

    print('Conver numpy massive to list of dicts...')

    new_edges = []
    nodes = [{'id': i, 'label': k, 'shape': 'image', 'size': 20, 'image': f'http://localhost:8000/{MEDIA_URL}icon.PNG'} for k, i in nodes.items()]
    for i in range(size):
        start_time = time.time()
        j = i
        while (j < size):
            if i != j:
                if edges[i][j] != 0:
                    new_edges.append({'from': i, 'to': j})
            j += 1
        if (i % 1000) == 0:
            print(f'Pass {i} row after {time.time() - start_time }')

    print('Done!')

    return nodes, new_edges



