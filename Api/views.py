import json
from flask import Blueprint, current_app, jsonify, request
from Bio import Entrez, Medline

from .model_analise import *
from .parsing import parse_record


analise = Blueprint('analise', __name__)


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
    print(len)
    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")
    records = [parse_record(record) for record in Medline.parse(handle) if not (record is None)]
    handle.close()

    current_app.logger.info(f'Start clear raw strings...')
    articles = create_clear_articles(records)

    current_app.logger.info(f'Start analise clear articles...')
    topics, props, embeddings = analise_articles(articles)
    records = return_results(records, topics, props)

    current_app.logger.info(f'Plot graphics for results...')
    current_app.logger.info(f'Plot clust graph')
    graph = return_clust_graph([rec['titl'] for rec in records], embeddings)

    current_app.logger.info(f'Plot heapmap')
    count_topics = len(set(topics))
    current_app.logger.info(f"Count topic = {count_topics}")
    n_clusters = 10
    if count_topics > 1:
        if n_clusters >= count_topics - 2:
            n_clusters = int(count_topics / 2)
        heapmap = return_heapmap(n_clusters=n_clusters)
        heapmap = json.loads(heapmap)
    else:
        heapmap = None

    current_app.logger.info(f'Plot heirarchy')
    try:
        heirarchy = return_heirarchy()
    except:
        heirarchy = None

    data = {
        'tematic_analise': records,
        'heapmap': heapmap,
        'clust_graph': json.loads(graph),
        'heirarchy': json.loads(heirarchy)
    }

    current_app.logger.info(f'Done analise!')

    return jsonify(data), 200