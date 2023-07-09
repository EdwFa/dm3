import requests


obj_color = {
    'disease': "#fdbbbb", # blue
    'drug': "#ECC58B",
    'gene': "#E2DB8C",
    'chemical': "#21c354",
    'species': "#A6EFDC",
    'mutation': "#B2DDEA",
    'cell_type': "#C6DEF5",
    'cell_line': "#A3B3D2",
    'DNA': "#C9B9E8",
    'RNA': "#D7DBE8",
}

def get_annotations(text, url="http://bern2.korea.ac.kr/plain"):
    print(len(text))
    if text == 'None AB':
        return None
    response = requests.post(url, json={'text': text})
    print(response.status_code)
    if response.status_code != 200:
        return None

    try:
        annotations = response.json()
    except Exception as e:
        print(e)
        print(response.content)
        return None
    return annotations


def get_pmid(pmid, url="http://bern2.korea.ac.kr/pubmed"):
    response = requests.get(f'{url}/{pmid}')
    if response.status_code != 200:
        print(f'Not exist {url}/{pmid}!')
        return None
    try:
        annotations = response.json()[0]
    except Exception as e:
        print(e)
        print(response.content)
        return None
    return annotations


def markup_text(text, annotations, **kwargs):
    markup_text = ''
    last_position = 0
    for annotation in annotations:
        print(annotation)
        start = annotation['span']['begin']
        end = annotation['span']['end']
        markup_str = f'<span style=\"color: {obj_color[annotation["obj"]]}\">{text[start:end]}<sub>{round(annotation["prob"], 2)}</sub></span>'
        markup_text = f'{markup_text}{text[last_position:start]}{markup_str}'

        last_position = end

    return markup_text
