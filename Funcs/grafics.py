import numpy
import time
import random

from affiliation_parser import parse_affil


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
        article_authors = article.auth.split('; ')
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
        authors = article.auth.split('; ')
        if len(authors) > author_on_article:
            authors = authors[:author_on_article]

        for author in authors:
            current_node = nodes[author]
            current_node.cluster_group = article.pl
            edges[current_node._id][current_node._id] += 1
            for a in authors:
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


def get_uniq_info_for_institutes(articles, author_on_article=10, count_of_rel=0):
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
    print('Get institutes...')
    for article in articles:
        article_authors = article.affl.split(';;; ')
        if len(article_authors) > author_on_article:
            article_authors = article_authors[:author_on_article]
        for aa in article_authors:
            authors.append(aa)


    authors = set(authors)
    size = len(authors)
    print('Count of institutes = ', size)
    print('Create nodes...')
    nodes = dict()
    for i, k in enumerate(authors):
        new_node = Node()
        new_node._id = i
        nodes[k] = new_node

    print('Getting info of edges...')
    edges = numpy.zeros((size, size), dtype='uint8')
    for article in articles:
        authors = article.affl.split(';;; ')
        if len(authors) > author_on_article:
            authors = authors[:author_on_article]
        used_affls = set()

        for author in authors:
            current_node = nodes[author]
            if not (author in used_affls):
                used_affls.add(author)
                current_node.cluster_group = article.pl
                edges[current_node._id][current_node._id] += 1
            for a in authors:
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
    for node in nodes:
        label_dict = parse_affil(node['label'])
        if label_dict['institution'] != '':
            node['label'] = label_dict['institution']
        else:
            print(label_dict)
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

def get_uniq_info_for_countries(articles, author_on_article=10, count_of_rel=0):
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
    print('Get counties...')
    countries = []
    for article in articles:
        article_authors = article.affl.split(';;; ')
        if len(article_authors) > author_on_article:
            article_authors = article_authors[:author_on_article]
        for aa in article_authors:
            authors.append(aa)
            if parse_affil(aa)['country'] == '':
                print(parse_affil(aa))
            countries.append(parse_affil(aa)['country'])


    countries = set(countries)
    size = len(countries)
    print('Count of countries = ', size)
    print('Create nodes...')
    nodes = dict()
    for i, k in enumerate(countries):
        new_node = Node()
        new_node._id = i
        nodes[k] = new_node


    print('Getting info of edges...')
    edges = numpy.zeros((size, size), dtype='uint8')
    for article in articles:
        authors = article.affl.split(';;; ')
        if len(authors) > author_on_article:
            authors = authors[:author_on_article]
        used_countries = set()

        for author in authors:
            author = parse_affil(author)['country']
            current_node = nodes[author]
            if not (author in used_countries):
                used_countries.add(author)
                edges[current_node._id][current_node._id] += 1
            for a in authors:
                a = parse_affil(a)['country']
                current_node2 = nodes[a]
                if current_node2._id == current_node._id:
                    continue
                edges[current_node._id][current_node2._id] += 1

    print(edges[:20, :20])
    print('Conver numpy massive to list of dicts...')

    new_edges = []
    random.seed(20)
    nodes = [{'id': node._id, 'label': k, 'x': random.uniform(-10., 10.), 'y': random.uniform(-10., 10.), 'weights': {'Documents': 1}} for k, node in nodes.items()]
    for node in nodes:
        if node['label'] == '':
            node['label'] = 'unknown'
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
    print('Get authors...')
    for article in articles:
        article_authors = article.auth.split('; ')
        for aa in article_authors:
            authors.append(aa)
        journals.append(article.jour)

    authors = set(authors)
    journals = set(journals)
    size = len(journals)
    print('Count of journals = ', size)
    print('Create nodes...')
    nodes = dict()
    for i, jorn in enumerate(journals):
        new_node = Node()
        new_node._id = i
        nodes[jorn] = new_node

    print('Getting info of journal`s edges...')
    edges = numpy.zeros((size, size), dtype='uint8')
    for author in authors:
        list_of_journals = []
        for article in articles:
            authors = article.auth.split('; ')
            if author in authors:
                list_of_journals.append(article.jour)
        list_of_journals = set(list_of_journals)
        for journal in list_of_journals:
            for j in list_of_journals:
                edges[nodes[journal]._id][nodes[j]._id] += 1

    print(edges[:20, :20])
    print('Convert numpy massive to list of dicts...')

    new_edges_1 = []
    random.seed(20)

    print('Convert journals...')
    nodes_1 = [{'id': node._id, 'label': k, 'x': random.uniform(-10., 10.), 'y': random.uniform(-10., 10.), 'weights': {'Documents': 1}} for k, node in nodes.items()]
    nodes_1.sort(key=sort_key)
    start_time = time.time()
    for i in range(size):
        j = i
        while (j < size):
            if i == j:
                nodes_1[i]['weights']['Documents'] = int(edges[i][i])
            else:
                if edges[i][j] > count_of_rel:
                    new_edges_1.append({'source_id': i, 'target_id': j, 'strength': int(edges[i][j])})

            j += 1
        if (i % 1000) == 0:
            print(f'Pass {i} row after {time.time() - start_time }')


    print('Done!')
    data = {
            'network': {
                'items': nodes_1,
                'links': new_edges_1,
            }
        }

    return data