from datetime import datetime, timezone, timedelta
from Bio import Entrez
from Bio import Medline
from celery import shared_task
import time
import json
import os
import math
import io
import http.client

import langchain
from langchain.chains import LLMChain
from langchain.schema import Document
from langchain.chat_models import ChatYandexGPT as YandexLLM
from langchain.callbacks import get_openai_callback
import grpc

import time
import jwt
import requests

from .serializers import *
from dm.settings import PARSER_EMAIL, MEDIA_URL, BERTTOPICAPI_URL, SUMMARISEAPI_URL, CHATAPI_URL, redis_cli, SERVICE_ACCOUNT_ID, KEY_ID, PRIVATE_KEY, IAM_TOKEN
from Funcs import *


def get_path_to_file(pk, file_name):
    path_to_file = os.path.join('datasets', str(pk), file_name)
    print(path_to_file)
    if not os.path.exists(path_to_file):
        raise Exception('Not found current_file')

    return path_to_file

def check_permission(user_permission, retmax):
    if user_permission.all_records is None:
        print('User has unlimited access!')
    else:
        allow_records = user_permission.all_records
        print(f'Allow analise records = {allow_records}, Len recieved data = {retmax}')
        if (retmax > allow_records):
            retmax = allow_records

    return retmax


@shared_task(bind=True)
def parse_records(self, query: str, count: int, new_task_id: int, retmax: int = 10000, email: str = None):
    # Парсим полученные записи
    print("Start parsing...")
    print(f"Count={count}, Task id={new_task_id}")
    new_task = TaskSearch.objects.get(id=new_task_id)
    print(f'Parsing = {retmax} records')

    if email:
        Entrez.email = email

    print(f'Email {email}')
    handle = Entrez.esearch(db="pubmed", sort='relevance', term=query, retmax=retmax)
    f = Entrez.read(handle)
    IdList = f['IdList']
    handle.close()

    i = 0
    records = []
    retries = 0

    while i < len(IdList):
        try:
            handle = Entrez.efetch(db="pubmed", id=IdList[i:i+1000], rettype="medline", retmode="text")
        except Exception as e:
            print(e)

            retries += 1
            if retries > 5:
                print('Max retries break')
                break
            time.sleep(1)
            print(f'Retry IdList from {i} to {i+1000}...')
            continue

        print(f'Parse IdList from {i} to {i+1000}')

        try:
            for record in Medline.parse(handle):
                data = parse_record(record)
                if data:
                    records.append(data)
                i += 1
        except Exception as e:
            print(f'Error on {i} step...')
            print(e)
            time.sleep(1)

        time.sleep(0.5)

        new_task.message = f'On {i}/{count} step...'
        new_task.save()

        handle.close()

    print(len(records))


    data = {
        'search_ncbi': ArticleSerializer(records, many=True).data
    }
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
def analise_records(self, pk, params, new_task_id):
    new_task = TaskAnalise.objects.get(id=new_task_id)
    new_task.message = 'Запущен процесс анализа, пожайлуста подождите...'
    new_task.save()

    params['current_task'] = new_task_id
    response = requests.post(f'{BERTTOPICAPI_URL}/analise', json=params, headers={'User-Agent': 'Mozilla/5.0'})
    print(response.status_code)
    if response.status_code != 200:
        raise Exception('Api is failed!')

    data = response.json()
    with open(get_path_to_file(pk, 'heapmap.json'), 'w') as f:
        json.dump(data['heapmap'], f)
    with open(get_path_to_file(pk, 'tematic_analise.json'), 'w') as f:
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
def summarise_text(self, records, email: str):
    response = requests.post(f'{SUMMARISEAPI_URL}/summarise', json={'IdList': records, 'email': email}, headers={'User-Agent': 'Mozilla/5.0'})
    print(response.status_code)
    if response.status_code != 200:
        raise Exception('Api is failed!')

    sentences = ' '.join(response.json()['summary'])

    return sentences


@shared_task(bind=True)
def summarise_emb(self, pk, text, email: str):
    text = ' '.join([s for s in text])
    response = requests.post(f'{SUMMARISEAPI_URL}/summarise_emb', json={'text': text, 'email': email}, headers={'User-Agent': 'Mozilla/5.0'})
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
def get_ddi_articles(self, query, new_task_id, **kwargs):
    url = f'https://www.ncbi.nlm.nih.gov/research/litsense-api/api/?query={query}&rerank=true'
    r = requests.get(url)
    print(r.status_code)
    if (r.status_code > 299 or r.status_code < 200):
        raise ConnectionError(f'Подключение к эмбедингам ncbi прошло неудачно с {r.status_code} ошибкой')

    records = json.loads(r.text)
    records_id = {record['pmid']: [round(record['score'], 2), record['section'], record['text']] for record in records
                  if current_record(record, **kwargs)}

    Entrez.email = kwargs.get('email', PARSER_EMAIL)
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
            record['use'] = True
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
def plot_graph_associations(self, IdList, pk, email, max_size=200):
    # Start plot graph on tematic analise

    if len(IdList) > max_size:
        IdList = IdList[:max_size - 1]

    Entrez.email = email
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

def get_token():
    try:
        # Получаем IAM-токен

        now = int(time.time())
        payload = {
            'aud': 'https://iam.api.cloud.yandex.net/iam/v1/tokens',
            'iss': SERVICE_ACCOUNT_ID,
            'iat': now,
            'exp': now + 360}

        # Формирование JWT
        encoded_token = jwt.encode(
            payload,
            PRIVATE_KEY,
            algorithm='PS256',
            headers={'kid': KEY_ID})

        url = 'https://iam.api.cloud.yandex.net/iam/v1/tokens'
        x = requests.post(url,
                          headers={'Content-Type': 'application/json'},
                          json={'jwt': encoded_token}).json()
        print(x)
        data = {
            'token': x['iamToken'],
            'time': x['expiresAt']
        }
        redis_cli.set(IAM_TOKEN, json.dumps(data))
        return True
    except Exception as e:
        print(e)
        return False

@shared_task(bind=True)
def translate(self, text):
    headers = {
        'Authorization': 'Api-Key AQVN02HWXvUTHSZwNaK_tsv3KpFmeBEQxDMJIFnJ',
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
    }
    body = {
        'targetLanguageCode': 'en',
        'format': "PLAIN_TEXT",
        'texts': [
            text
        ],
        'speller': True
    }
    response = requests.post('https://translate.api.cloud.yandex.net/translate/v2/translate', headers=headers,
                             json=body)

    print(response.status_code)

    return response.json()

def create_doc(record, i):
    return Document(page_content=f"{record}\n", metadata={'page': f"{i}"})

@shared_task(bind=True)
def send_message(self, **kwargs):
    records = [create_doc(record, i) for i, record in enumerate(kwargs['articles'])]
    len_records = 0
    for record in records:
        len_records += len(record.page_content.split(' '))
    print(len_records)
    if not redis_cli.exists(IAM_TOKEN):
        if not get_token():
            raise Exception('Token didnt got')

    token_data = json.loads(redis_cli.get(IAM_TOKEN))
    token = token_data['token']
    print(token_data)
    try:
        token_time = datetime.strptime(token_data['time'], '%m/%d/%Y:%H-%M')
    except:
        token_time = datetime.strptime(token_data['time'].split('.')[0], '%Y-%m-%dT%H:%M:%S')
    timed = token_time - datetime.now()
    if timed.days < 0:
        timed = -1 * timed.days * 24 - timed.seconds / 3600
    else:
        timed = timed.days * 24 + timed.seconds / 3600
    print(token_time, datetime.now(), timed)

    if timed > 1:
        get_token()

    instructions = """Представь себе, что ты сотрудник Yandex Cloud. Твоя задача - вежливо и по мере своих сил отвечать на все вопросы собеседника."""

    llm = YandexLLM(iam_token=token,
                    instruction_text=instructions,
                    timeout=30)
    # Промпт для обработки документов
    document_prompt = langchain.prompts.PromptTemplate(
        input_variables=["page_content"],
        template="{page_content}"
    )

    # Промпт для языковой модели
    document_variable_name = "context"
    stuff_prompt_override = """
        Пожалуйста, посмотри на текст ниже и ответь на вопрос, используя информацию из этого текста или, если нет информации, своими словами ответь что знаешь.
        Текст:
        -----
        {context}
        -----
        Вопрос:
        {query}
    """
    prompt = langchain.prompts.PromptTemplate(
        template=stuff_prompt_override,
        input_variables=["context", "query"]
    )

    # Создаём цепочку
    llm_chain = langchain.chains.LLMChain(llm=llm,
                                          prompt=prompt)

    chain = langchain.chains.combine_documents.stuff.StuffDocumentsChain(
        llm_chain=llm_chain,
        document_prompt=document_prompt,
        document_variable_name=document_variable_name,
    )
    try:
        message = chain.run(input_documents=records, query=kwargs['message'])
        print(message)
    except grpc.RpcError as e:
        print(e)
        message = "Слишком большой запрос по кол-ву документов, пожайлуста уменьшите их кол-во и попробуйте снова"

    return message

