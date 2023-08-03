import json
from flask import Blueprint, current_app, jsonify, request

from .YandexChatGPT import model, getResponse


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
    if not ('message' in data):
        return jsonify({'status': 'Error', 'message': 'message not received'}), 500

    response = getResponse(data['message'])

    data = {
        'message': response
    }

    current_app.logger.info(f'Done response!')

    return jsonify(data), 200