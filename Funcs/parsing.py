from datetime import date
from Bio import Entrez

from main.models import Article


# Запрос для начального запуска как пример
default_query = "{0} and {1}:{2}[dp] and Clinical Trial[pt]".format('covid-19', '2000/01/01',
                                                                    date.today().strftime("%Y/%m/%d"))

# Создание запроса на основе полученных ключей
def create_query(**filters):
    # Преобразуем наши ключи в запрос pubmed
    print([i for i in filters.items()])

    if len([i for i in filters]) == 0:
        return default_query

    # Содаем запрос с текстом запроса
    query = ""
    text = filters.get('search_field', None)
    if text:
        query = f'{text}'

    # Добавляем даты
    start_date = filters.get('dateStart', "1900-01-01")
    end_date = filters.get('dateStop', date.today().strftime("%Y-%m-%d"))
    query = f"{query} and {start_date.replace('-', '/')}:{end_date.replace('-', '/')}[dp]"

    # Добавляем гендеры
    genders = filters.get('Gender', [])
    if len(genders) != 0:
        query_genders = " OR ".join([f'{i}[mh]' for i in genders])
        query = f"{query} AND ({query_genders})"

    # Добавляем типы документов
    types = filters.get('Type', [])
    if len(types) != 0:
        query_types = " OR ".join([f'{i}[pt]' for i in types])
        query = f"{query} AND ({query_types})"

    # Добавляем Возрастные группы
    olds = filters.get('Old', [])
    if len(olds) != 0:
        query_olds = " OR ".join(olds)
        query = f"{query} AND ({query_olds})"

    print(query)

    return query


def get_records(query: str):
    # Получаем наши записи
    Entrez.email = "e.p@d_health.pro"  # Говорю NCBI кто я есть
    handle = Entrez.esearch(db="pubmed", sort='relevance', term=query)
    record = Entrez.read(handle)

    handle.close()

    str_query = record['QueryTranslation']

    try:
        translation_stack = record['TranslationStack']
    except:
        translation_stack = 'null/null/null'

    records_count = int(record['Count'])

    print(f"All records: {records_count}")

    return str_query, translation_stack, records_count


# Парсинг статьи
def parse_record(record):
    # Парсим полученный словарь записи
    if not ('PMID' in record): # Проверяем на наличие pmid если его нет не сохраняем данные
        return

    data = Article()
    data.uid = record['PMID']
    data.url = 'https://pubmed.ncbi.nlm.nih.gov/' + record['PMID'] + '/'
    if 'TI' in record:
        data.titl = record['TI']
    else:
        data.titl = "None Tilte"

    if 'AID' in record:
        data.aid = ' '.join([''.join(i.split(' ')[:-1]) for i in record['AID']])
        data.urlaid = ' '.join([f'https://sci-hub.do/{aid}' for aid in data.aid.split(' ')])
    else:
        data.aid = "Без AID"
        data.urlaid = "Без url AID"

    if 'JT' in record:
        data.jour = record['JT']
    else:
        data.jour = "без JT"
    if 'DP' in record:
        data.pdat = record['DP']
    else:
        data.pdat = "без DP"
    if 'FAU' in record:
        data.auth = '; '.join(record['FAU'])
    else:
        data.auth = "не указаны"

    if 'AD' in record:
        data.affl = ';;; '.join(record['AD'])
    else:
        data.affl = "None AD"

    if 'PT' in record:
        data.pt = '; '.join(record['PT'])
    else:
        data.pt = "без PT"

    if 'PL' in record:
        data.pl = record['PL']
    else:
        data.pl = "без PL"

    if 'MH' in record:
        data.mesh = ';;; '.join(record['MH'])
    else:
        data.mesh = "без MH"

    if 'AB' in record:
        data.tiab = record['AB']
    else:
        data.tiab = "None AB"

    if 'OT' in record:
        data.ot = ';;; '.join(record['OT'])
    else:
        data.ot = "без OT"

    if 'RN' in record:
        data.rn = ';;; '.join(record['RN'])
    else:
        data.rn = "без RN"

    return data