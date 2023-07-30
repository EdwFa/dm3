import json
from flask import Blueprint, current_app, jsonify, request
from flask import send_file
import nltk
from Bio import Entrez, Medline

from .summarise import *
from .parsing import parse_record


sum = Blueprint('sum', __name__)


@sum.route('/check', methods=['GET'])
def check_status():
    return jsonify({'status': 'Work'}), 200


@sum.route('/summarise', methods=['POST'])
def summarise_records():
    """
    -> Receive pmids from main server
    -> Get all article from pmid
    -> Analise articles in BertTopic
    -> Use learned model for print graphics
    -> """
    current_app.logger.info(f'Start summarise text...')

    data = request.json
    if not ('IdList' in data):
        return jsonify({'status': 'Error', 'message': 'Text not found'}), 500
    IdList = data['IdList']
    handle = Entrez.efetch(db="pubmed", id=IdList, rettype="medline", retmode="text")
    records = [parse_record(record) for record in Medline.parse(handle) if not (record is None)]
    raw_text = ' '.join([rec['titl'] + rec['tiab'] for rec in records])
    handle.close()

    current_app.logger.info(f'Делим документ на предложения...')
    sentences = nltk.sent_tokenize(raw_text)
    print("Num sentences:", len(sentences))

    current_app.logger.info(f'Формируем эмбединги по предложениям...')
    embeddings = model.encode(sentences, convert_to_tensor=True)

    current_app.logger.info(f'Считаем косинусную близость для пар предложений...')
    cos_scores = util.cos_sim(embeddings, embeddings).cpu().numpy()

    current_app.logger.info(f'Центральные метрики предложений...')
    centrality_scores = degree_centrality_scores(cos_scores, threshold=None)

    current_app.logger.info(f'Сортируем предложения с наивысшими далами в начало списка для последующих срезов...')
    most_central_sentence_indices = np.argsort(-centrality_scores)

    current_app.logger.info(f'Саммари из 10 предложений...')
    most_central_sentences = []
    print("\n\nSummary:")
    for idx in most_central_sentence_indices[0:10]:
        print(sentences[idx].strip())
        most_central_sentences.append(sentences[idx].strip())

    data = {
        'status': 'Done!',
        'summary': most_central_sentences
    }

    current_app.logger.info(f'Done analise!')

    return jsonify(data), 200

@sum.route('/summarise_emb', methods=['POST'])
def summarise_emb():
    """
    -> Receive pmids from main server
    -> Get all article from pmid
    -> Analise articles in BertTopic
    -> Use learned model for print graphics
    -> """
    current_app.logger.info(f'Start summarise text...')

    data = request.json
    if not ('text' in data):
        return jsonify({'status': 'Error', 'message': 'Text not found'}), 500
    raw_text = data['text']

    current_app.logger.info(f'Делим документ на предложения...')
    sentences = nltk.sent_tokenize(raw_text)
    print("Num sentences:", len(sentences))

    current_app.logger.info(f'Формируем эмбединги по предложениям...')
    embeddings = model.encode(sentences, convert_to_tensor=True)

    current_app.logger.info(f'Считаем косинусную близость для пар предложений...')
    cos_scores = util.cos_sim(embeddings, embeddings).numpy()

    current_app.logger.info(f'Центральные метрики предложений...')
    centrality_scores = degree_centrality_scores(cos_scores, threshold=None)

    current_app.logger.info(f'Сортируем предложения с наивысшими далами в начало списка для последующих срезов...')
    most_central_sentence_indices = np.argsort(-centrality_scores)

    current_app.logger.info(f'Саммари из 10 предложений...')
    most_central_sentences = []
    print("\n\nSummary:")
    for idx in most_central_sentence_indices[0:10]:
        print(sentences[idx].strip())
        most_central_sentences.append(sentences[idx].strip())

    data = {
        'status': 'Done!',
        'summary': most_central_sentences
    }

    current_app.logger.info(f'Done analise!')

    return jsonify(data), 200
