import numpy
import time
import random


def sort_key(d):
    return d['id']


class NodeAuthor:
    _id = 0
    cluster_group = 'null'

    def __str__(self):
        return f'{self._id} - {self.cluster_group}'


def get_uniq_info_for_graph(articles, author_on_article=10, count_of_rel=0):
    """
    return json file formated
    data = {
        'networks': {
            "items": [
                {
                    "id": pk :int,
                    "label": author name :str,
                    "x": x coord :float,
                    "y": y coord :float,
                    "weights": {
                        "Documents": count of apperance: int,
                    },
                },
                ...
            ]
            "links": [
                {
                    "source_id": id item first: pk int,
                    "target_id": id item second: pk int,
                    "strength": count of pub together: int > 0,
                },
                ...
            ]
        }
    }
    """
    authors = []
    print('Get authors...')
    for article in articles:
        article_authors = article['auth'].split('; ')
        if len(article_authors) > author_on_article:
            article_authors = article_authors[:author_on_article]
        for aa in article_authors:
            authors.append(aa)

    authors = set(authors)
    size = len(authors)
    print('Count of authors = ', size)
    print('Create nodes...')
    nodes = dict()
    for i, k in enumerate(authors):
        new_node = NodeAuthor()
        new_node._id = i
        nodes[k] = new_node

    print('Getting info of edges...')
    edges = numpy.zeros((size, size), dtype='uint8')
    for article in articles:
        authors = article['auth'].split('; ')
        if len(authors) > author_on_article:
            authors = authors[:author_on_article]
        for author in authors:
            current_node = nodes[author]
            current_node.cluster_group = article['pl']
            for a in authors:
                edges[current_node._id][nodes[a]._id] += 1

    print(edges[:20, :20])
    print('Conver numpy massive to list of dicts...')

    new_edges = []
    random.seed(20)
    uniq_contries = {k: v + 1 for v, k in enumerate(set([i.cluster_group for i in nodes.values()]))}
    nodes = [{'id': node._id, 'label': k, 'x': random.uniform(-10., 10.), 'y': random.uniform(-10., 10.), 'cluster': uniq_contries[node.cluster_group], 'weights': {'Documents': 1}} for k, node in nodes.items()]
    nodes.sort(key=sort_key)
    start_time = time.time()
    for i in range(size):
        j = i
        while (j < size):
            if i == j:
                nodes[i]['weights']['Documents'] = int(edges[i][i])
            else:
                if edges[i][j] > count_of_rel:
                    new_edges.append({'source_id': i, 'target_id': j, 'strength': int(edges[i][j])})

            j += 1
        if (i % 1000) == 0:
            print(f'Pass {i} row after {time.time() - start_time }')

    print('Done!')
    data = {
        'network': {
            'items': nodes,
            'links': new_edges,
        }
    }
    return data