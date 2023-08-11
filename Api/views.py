import json
import os

from flask import Blueprint, current_app, jsonify, request
from flask import send_file
from Bio import Entrez, Medline
import time
import redis

from .model_analise import *
from .parsing import parse_record


analise = Blueprint('analise', __name__)

redis_cli = redis.Redis(host=os.getenv('REDIS_HOST'), port=os.getenv('REDIS_PORT'), decode_responses=True)

@analise.route('/check', methods=['GET'])
def check_status():
    return jsonify({'status': 'Work'}), 200


@analise.route('/analise', methods=['POST'])
def analise_records():
    """
    -> Receive pmids from main server
    -> Get all article from pmid
    -> Analise articles in BertTopic
    -> Use learned model for print graphics
    -> """
    current_app.logger.info(f'Start analise records...')

    current_app.logger.info(f'Receive data(articlesId)...')
    data = request.json
    if not ('articles' in data):
        return jsonify({'status': 'Error', 'message': 'No one article id'}), 500
    IdList = data['articles']

    current_app.logger.info(f'Get limit on analise package ({len(IdList)})...')
    email = data.get('email', None)
    if email:
        Entrez.email = email
    filters = data['filters']

    current_task = data['current_task']

    #redis_cli.get(f'task:{current_task}')

    i = 0
    records = []
    retries = 0

    while i < len(IdList):
        redis_cli.set(f'task:{current_task}', 'Передаем данные на анализ...')
        try:
            handle = Entrez.efetch(db="pubmed", id=IdList[i:i + 1000], rettype="medline", retmode="text")
        except Exception as e:
            print(e)

            retries += 1
            if retries > 5:
                print('Max retries break')
                break
            time.sleep(1)
            print(f'Retry IdList from {i} to {i + 1000}...')
            continue

        print(f'Parse IdList from {i} to {i + 1000}')

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

        handle.close()

    # handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")
    # records = [parse_record(record) for record in Medline.parse(handle) if not (record is None)]
    # handle.close()

    current_app.logger.info(f'Start clear raw strings...')
    redis_cli.set(f'task:{current_task}', ' Токенизация и очистка от стоп слов...')
    articles = create_clear_articles(records)

    current_app.logger.info(f'Start analise clear articles...')
    redis_cli.set(f'task:{current_task}', 'Построение построение эмбедингов...')
    current_app.logger.info(f'Bild Bert...')
    current_app.logger.info(filters)
    topic_model = bild_bert(**filters)
    topics, props, embeddings = analise_articles(topic_model, articles)
    records = return_results(records, topics, props)

    current_app.logger.info(f'Plot graphics for results...')
    current_app.logger.info(f'Plot clust graph')
    redis_cli.set(f'task:{current_task}', 'Cжатие по UMAP...')
    graph = return_clust_graph(topic_model, [rec['titl'] for rec in records], embeddings, **filters)

    current_app.logger.info(f'Plot heapmap')
    redis_cli.set(f'task:{current_task}', 'Отрисовка графиков...')
    count_topics = len(set(topics))
    current_app.logger.info(f"Count topic = {count_topics}")
    n_clusters = int(filters.get('n_clusters', 10))
    if count_topics > 1:
        if n_clusters >= count_topics - 2:
            n_clusters = int(count_topics / 2)
        heapmap = return_heapmap(topic_model, top_n_topics=int(filters.get('top_n_topics', 40)), n_clusters=n_clusters)
        heapmap = json.loads(heapmap)
    else:
        heapmap = None

    current_app.logger.info(f'Plot heirarchy')
    try:
        heirarchy = json.loads(return_heirarchy(topic_model))
    except:
        heirarchy = None

    current_app.logger.info(f'Plot DTM')
    try:
        DTM = json.loads(return_DTM(topic_model, [rec['titl'] for rec in records], [int(record['pdat'].split(' ')[0]) for record in records]))
    except:
        DTM = None

    data = {
        'tematic_analise': records,
        'heapmap': heapmap,
        'clust_graph': json.loads(graph),
        'heirarchy': heirarchy,
        'DTM': DTM,
        'embeddings': embeddings.tolist(),
        'topics': return_topic_label(topic_model)
    }

    redis_cli.set(f'task:{current_task}', 'Передача назад...')
    current_app.logger.info(f'Done analise!')

    return jsonify(data), 200

@analise.route('/download_vectors', methods=['GET'])
def downloadVerctors():
    data = request.json
    if 'user' not in data:
        return jsonify({'status': 'Error', 'message': 'No found user!'}), 500
    user_id = data['user'].id
    path = f"/{user_id}/vectors_OR.tsv"
    return send_file(path, as_attachment=True)

@analise.route('/download_vectors', methods=['GET'])
def downloadMetaData():
    data = request.json
    if 'user' not in data:
        return jsonify({'status': 'Error', 'message': 'No found user!'}), 500
    user_id = data['user'].id
    path = f"/{user_id}/metadata_OR.tsv"
    return send_file(path, as_attachment=True)