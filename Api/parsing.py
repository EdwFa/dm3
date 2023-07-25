class Article:
    uid = 0
    aid = ""
    titl = ""
    mesh = ""
    majr = ""
    subh = ""
    auth = ""
    jour = ""
    affl = ""
    pdat = ""
    tiab = ""
    ptyp = ""
    url = ""
    urlaid = ""
    pt = ""
    pl = ""


# Парсинг статьи
def parse_record(record):
    # Парсим полученный словарь записи
    if not ('PMID' in record): # Проверяем на наличие pmid если его нет не сохраняем данные
        return None

    data = Article()
    data.uid = int(record['PMID'])
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

    return data.__dict__