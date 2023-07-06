import numpy
import time


def get_uniq_info_for_graph(articles, author_on_article=10, count_of_rel=0):
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
                if edges[i][j] > count_of_rel:
                    new_edges.append({'from': i, 'to': j})
            j += 1
        if (i % 1000) == 0:
            print(f'Pass {i} row after {time.time() - start_time }')

    print('Done!')

    return nodes, new_edges