import numpy
import time
import random


def sort_key(d):
    return d['id']


class Node:
    _id = 0
    cluster_group = 'null'

    def __str__(self):
        return f'{self._id} - {self.cluster_group}'


def get_uniq_info_for_authors(articles, author_on_article=10, count_of_rel=0):
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
    print('Get authors, journals and countries...')
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
        new_node = Node()
        new_node._id = i
        nodes[k] = new_node

    print('Getting info of edges...')
    edges = numpy.zeros((size, size), dtype='uint8')
    for article in articles:
        authors = article['auth'].split('; ')
        if len(authors) > author_on_article:
            authors = authors[:author_on_article]

        if len(authors) > 4:
            main_authors = authors[:3]
        else:
            main_authors = authors

        for author in authors:
            current_node = nodes[author]
            current_node.cluster_group = article['pl']
            edges[current_node._id][current_node._id] += 1
            for a in main_authors:
                current_node2 = nodes[a]
                if current_node2._id == current_node._id:
                    continue
                edges[current_node._id][current_node2._id] += 1

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


def get_uniq_info_for_other(articles, count_of_rel=0):
    """
    return json file formated for journals
    data = {
        'networks': {
            "items": [
                {
                    "id": pk :int,
                    "label": journal name :str,
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
                    "strength": count of pub together there same author: int > 0,
                },
                ...
            ]
        }
    }
    """
    authors = []
    journals = []
    countries = []
    print('Get authors...')
    for article in articles:
        article_authors = article['auth'].split('; ')
        for aa in article_authors:
            authors.append(aa)
        journals.append(article['jour'])
        countries.append(article['pl'])

    authors = set(authors)
    journals = set(journals)
    countries = set(countries)
    size_1 = len(journals)
    size_2 = len(countries)
    print('Count of journals = ', size_1)
    print('Count of countries = ', size_2)
    print('Create nodes...')
    nodes_1 = dict()
    for i, k in enumerate(journals):
        new_node = Node()
        new_node._id = i
        nodes_1[k] = new_node
    nodes_2 = dict()
    for i, k in enumerate(countries):
        new_node = Node()
        new_node._id = i
        nodes_2[k] = new_node

    print('Getting info of journal`s edges...')
    edges_1 = numpy.zeros((size_1, size_1), dtype='uint8')
    edges_2 = numpy.zeros((size_2, size_2), dtype='uint8')
    for author in authors:
        list_of_journals = []
        list_of_countries = []
        for article in articles:
            authors = article['auth'].split('; ')
            if author in authors:
                list_of_journals.append(article['jour'])
                list_of_countries.append(article['pl'])
        list_of_journals = set(list_of_journals)
        list_of_countries = set(list_of_countries)
        for journal in list_of_journals:
            for j in list_of_journals:
                edges_1[nodes_1[journal]._id][nodes_1[j]._id] += 1
        for country in list_of_countries:
            for c in list_of_countries:
                edges_2[nodes_2[country]._id][nodes_2[c]._id] += 1

    print(edges_1[:20, :20])
    print(edges_2[:20, :20])
    print('Convert numpy massive to list of dicts...')

    new_edges_1 = []
    new_edges_2 = []
    random.seed(20)

    print('Convert journals...')
    nodes_1 = [{'id': node._id, 'label': k, 'x': random.uniform(-10., 10.), 'y': random.uniform(-10., 10.), 'weights': {'Documents': 1}} for k, node in nodes_1.items()]
    nodes_1.sort(key=sort_key)
    start_time = time.time()
    for i in range(size_1):
        j = i
        while (j < size_1):
            if i == j:
                nodes_1[i]['weights']['Documents'] = int(edges_1[i][i])
            else:
                if edges_1[i][j] > count_of_rel:
                    new_edges_1.append({'source_id': i, 'target_id': j, 'strength': int(edges_1[i][j])})

            j += 1
        if (i % 1000) == 0:
            print(f'Pass {i} row after {time.time() - start_time }')

    print('Convert countries...')
    nodes_2 = [{'id': node._id, 'label': k, 'x': random.uniform(-10., 10.), 'y': random.uniform(-10., 10.), 'weights': {'Documents': 1}} for k, node in nodes_2.items()]
    nodes_2.sort(key=sort_key)
    start_time = time.time()
    for i in range(size_2):
        j = i
        while (j < size_2):
            if i == j:
                nodes_2[i]['weights']['Documents'] = int(edges_2[i][i])
            else:
                if edges_2[i][j] > count_of_rel:
                    new_edges_2.append({'source_id': i, 'target_id': j, 'strength': int(edges_2[i][j])})

            j += 1
        if (i % 1000) == 0:
            print(f'Pass {i} row after {time.time() - start_time}')


    print('Done!')
    data = [
        {
            'network': {
                'items': nodes_1,
                'links': new_edges_1,
            }
        },
        {
            'network': {
                'items': nodes_2,
                'links': new_edges_2,
            }
        }
    ]
    return data