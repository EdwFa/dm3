from nltk.tokenize import word_tokenize
import numpy
from umap import UMAP
from bertopic import BERTopic
import plotly
from sentence_transformers import SentenceTransformer


EMBEDDING_MODEL = 'multi-qa-MiniLM-L6-cos-v1'
_min_TOPIC_SIZE = 10
_top_N_WORDS = 15
_n_gram_range = (1, 3)

change_simpols = ['"', "'", '[', ']', '\\', '/']
stop_words = ['i', 'me', 'my', 'myself', 'we', 'our', 'ours', 'ourselves',
              'you', "you're", "you've", "you'll", "you'd", 'your', 'yours', 'yourself', 'yourselves',
              'he', 'him', 'his', 'himself', 'she', "she's", 'her', 'hers', 'herself',
              'it', "it's", 'its', 'itself',
              'they', 'them', 'their', 'theirs', 'themselves',
              'what', 'which', 'who', 'whom',
              'this', 'that', "that'll", 'these', 'those',
              'am', 'is', 'are', 'was', 'were', 'be', 'been', 'being',
              'have', 'has', 'had', 'having',
              'do', 'does', 'did', 'doing',
              'a', 'an', 'the', 'and', 'but', 'if',
              'or', 'because', 'as',
              'until', 'while', 'of', 'at', 'by', 'for', 'with', 'about',
              # 'against', 'between', 'into', 'through', 'during', 'before', 'after', 'above', 'below',
              'to', 'from', 'up', 'down', 'in', 'out', 'on', 'off', 'over', 'under', 'again', 'further',
              'then', 'once', 'here', 'there', 'when', 'where', 'why', 'how', 'all', 'any', 'both', 'each', 'few', 'more', 'most', 'other', 'some', 'such',
              'no', 'nor', 'not', 'only', 'own', 'same', 'so', 'than', 'too', 'very', 's', 't', 'can', 'will', 'just', 'don', "don't", 'should', "should've", 'now',
              'd', 'll', 'm', 'o', 're', 've', 'y', 'ain', 'aren', "aren't", 'couldn', "couldn't", 'didn', "didn't", 'doesn', "doesn't", 'hadn', "hadn't", 'hasn', "hasn't",
              'haven', "haven't", 'isn', "isn't", 'ma', 'mightn', "mightn't", 'mustn', "mustn't", 'needn', "needn't", 'shan', "shan't", 'shouldn', "shouldn't",
              'wasn', "wasn't", 'weren', "weren't", 'won', "won't", 'wouldn', "wouldn't"]


def bild_bert(**paramets):
    n_gram_range = (int(paramets.get('rangeMin', 1)), int(paramets.get('rangeMax', 3)))
    min_TOPIC_SIZE = int(paramets.get('min_TOPIC_SIZE', _min_TOPIC_SIZE))
    top_N_WORDS = int(paramets.get('top_N_WORDS', _top_N_WORDS))
    topic_model = BERTopic(embedding_model=EMBEDDING_MODEL, n_gram_range=n_gram_range, min_topic_size=min_TOPIC_SIZE,
                           top_n_words=top_N_WORDS, verbose=True)
    return topic_model

def clear_text(text):
    for change_simbol in change_simpols:
        text = text.replace(change_simbol, ' ')
    tokens = word_tokenize(text)
    text = ' '.join([word for word in tokens if not word.lower() in stop_words])
    return text

def create_clear_articles(records):
    clear_articles = numpy.array(['no one' for i in records], dtype=object)
    for i, record in enumerate(records):
        print("RECORD ===>")
        print(record)
        clear_articles[i] = clear_text(record['titl'] + " " + record['tiab'])
    return clear_articles

def analise_articles(topic_model, articles):
    sentence_model = SentenceTransformer(EMBEDDING_MODEL)
    embeddings = sentence_model.encode(articles, show_progress_bar=True)
    # Create topic model
    topics, probs = topic_model.fit_transform(articles, embeddings)
    return topics, probs, embeddings

def return_results(records, topics, props):
    for rec, topic, prop in zip(records, topics, props):
        rec['topic'] = topic
        rec['prop'] = round(prop, 2)
    return records

def return_clust_graph(topic_model, docs_title, embeddings, **params):
    n_neighbors = int(params.get('n_neighbors', 10))
    n_components = int(params.get('n_components', 2))
    min_dist = float(params.get('min_dist', 0.0))
    metric = params.get('metric', 'cosine')
    reduced_embeddings = UMAP(n_neighbors=n_neighbors, n_components=n_components, min_dist=min_dist, metric=metric).fit_transform(embeddings)
    fig = topic_model.visualize_documents(docs_title, embeddings=reduced_embeddings, width=1000, height=1000)
    graphJSON = plotly.io.to_json(fig, pretty=True)
    return graphJSON

def return_heapmap(topic_model, top_n_topics: int, n_clusters: int):
    fig = topic_model.visualize_heatmap(top_n_topics=top_n_topics, n_clusters=n_clusters, width=1000, height=1000, custom_labels=True)
    graphJSON = plotly.io.to_json(fig, pretty=True)
    return graphJSON

def return_heirarchy(topic_model):
    fig = topic_model.visualize_hierarchy()
    graphJSON = plotly.io.to_json(fig, pretty=True)
    return graphJSON

def return_DTM(topic_model, docs_title, dates):
    topics_over_time = topic_model.topics_over_time(docs_title, dates)
    fig = topic_model.visualize_topics_over_time(topics_over_time)
    graphJSON = plotly.io.to_json(fig, pretty=True)
    return graphJSON