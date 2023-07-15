import os

from flask import Flask
import nltk
import logging

from Api import analise, check_device, Entrez


app = Flask(__name__, static_url_path='')
app.register_blueprint(analise)
app.logger.name = 'BertTopicApi'
Entrez.email = os.getenv('PARSER_EMAIL')  # Говорю NCBI кто я есть
address_host = os.getenv('HOST')
address_port = os.getenv('PORT')


nltk.download('stopwords')
nltk.download('punkt')


if __name__ == '__main__':
    logging.basicConfig(filename='INFO.log', level=logging.DEBUG)
    app.logger.info('-----------------')
    app.logger.info('Starting...')
    app.logger.info(f'GPU in system = {check_device()}')
    app.run(debug=True, host=address_host, port=address_port)
