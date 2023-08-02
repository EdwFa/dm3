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
    data = Article()

    if 'TI' in record:
        data.titl = record['TI']
    else:
        data.titl = "None Title"

    if 'AB' in record:
        data.tiab = record['AB']
    else:
        data.tiab = "None AB"

    return data