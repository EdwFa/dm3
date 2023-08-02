from django.db import models
from django.contrib.auth import get_user_model

from datetime import datetime

User = get_user_model()


task_statuses = (
    (0, 'processing'),
    (1, 'done'),
    (2, 'error'),
)

task_analise = (
    (0, 'tematic'),
    (1, 'markup'),
    (2, 'graph'),
)


# Create your models here.
class Article(models.Model):
    uid = models.IntegerField('PMID', primary_key=True)
    aid = models.CharField('AID', max_length=64)
    titl = models.CharField('TITL', max_length=250)
    mesh = models.TextField('MESH')
    majr = models.TextField('MAJR')
    subh = models.TextField('SUBH')
    auth = models.TextField('AUTH')
    jour = models.TextField('JOUR')
    affl = models.TextField('AFFL')
    pdat = models.CharField('PDAT', max_length=32)
    tiab = models.TextField('TIAB')
    ptyp = models.TextField('PTYP')
    url = models.CharField('URL_LINK', max_length=250)
    urlaid = models.CharField('URL_AID', max_length=250)
    pt = models.TextField('PT')
    pl = models.TextField('PL')

    def __str__(self):
        return str(self.uid)

'''
Структура и поля публикаций pubmed
ALL, All Fields, All terms from all searchable fields
UID, UID, Unique number assigned to publication
FILT, Filter, Limits the records
TITL, Title, Words in title of publication
WORD, Text Word, Free text associated with publication
MESH, MeSH Terms, Medical Subject Headings assigned to publication
MAJR, MeSH Major Topic, MeSH terms of major importance to publication
AUTH, Author, Author(s) of publication
JOUR, Journal, Journal abbreviation of publication
AFFL, Affiliation, Author's institutional affiliation and address
ECNO, EC/RN Number, EC number for enzyme or CAS registry number
SUBS, Supplementary Concept, CAS chemical name or MEDLINE Substance Name
PDAT, Date - Publication, Date of publication
EDAT, Date - Entrez, Date publication first accessible through Entrez
VOL, Volume, Volume number of publication
PAGE, Pagination, Page number(s) of publication
PTYP, Publication Type, Type of publication (e.g., review)
LANG, Language, Language of publication
ISS, Issue, Issue number of publication
SUBH, MeSH Subheading, Additional specificity for MeSH term
SI, Secondary Source ID, Cross-reference from publication to other databases
MHDA, Date - MeSH, Date publication was indexed with MeSH terms
TIAB, Title/Abstract, Free text associated with Abstract/Title
OTRM, Other Term, Other terms associated with publication
INVR, Investigator, Investigator
COLN, Author - Corporate, Corporate Author of publication
CNTY, Place of Publication, Country of publication
PAPX, Pharmacological Action, MeSH pharmacological action pre-explosions
GRNT, Grant Number, NIH Grant Numbers
MDAT, Date - Modification, Date of last modification
CDAT, Date - Completion, Date of completion
PID, Publisher ID, Publisher ID
FAUT, Author - First, First Author of publication
FULL, Author - Full, Full Author Name(s) of publication
FINV, Investigator - Full, Full name of investigator
TT, Transliterated Title, Words in transliterated title of publication
LAUT, Author - Last, Last Author of publication
PPDT, Print Publication Date, Date of print publication
EPDT, Electronic Publication Date, Date of Electronic publication
LID, Location ID, ELocation ID
CRDT, Date - Create, Date publication first accessible through Entrez
BOOK, Book, ID of the book that contains the document
ED, Editor, Section's Editor
ISBN, ISBN, ISBN
PUBN, Publisher, Publisher's name
AUCL, Author Cluster ID, Author Cluster ID
EID, Extended PMID, Extended PMID
DSO, DSO, Additional text from the summary
AUID, Author - Identifier, Author Identifier
PS, Subject - Personal Name, Personal Name as Subject
COIS, Conflict of Interest Statements, Conflict of Interest Statements
'''


class TaskSearch(models.Model):
    task_id = models.CharField(max_length=150, null=True)
    user = models.ForeignKey(User, on_delete=models.DO_NOTHING, null=True, related_name='task_search')
    query = models.TextField()
    full_query = models.TextField()
    translation_stack = models.TextField()
    start_date = models.DateTimeField(auto_now_add=True)
    end_date = models.DateTimeField(null=True, blank=True)
    status = models.IntegerField(choices=task_statuses, default=0)
    count = models.IntegerField()
    message = models.TextField(default='')

    def __str__(self):
        return f'{self.user.email}: {self.query} ---> ({task_statuses[self.status][1]})'

    def get_user(self):
        return self.user.email

    def get_status(self):
        return task_statuses[self.status][1]

    def time_delta(self):
        if self.end_date is None:
            return None
        delta_time = (self.end_date - self.start_date).seconds
        return delta_time


class TaskAnalise(models.Model):
    task_id = models.CharField(max_length=150, null=True)
    user = models.ForeignKey(User, on_delete=models.DO_NOTHING, null=True, related_name='task_analise')
    start_date = models.DateTimeField(auto_now_add=True)
    end_date = models.DateTimeField(null=True, blank=True)
    status = models.IntegerField(choices=task_statuses, default=0)
    type_analise = models.IntegerField(choices=task_analise, default=0)
    message = models.TextField(default='')

    def __str__(self):
        return f'{self.user.email}: {task_analise[self.type_analise][1]} ---> ({task_statuses[self.status][1]})'

    def get_user(self):
        return self.user.email

    def get_type_analise(self):
        if self.type_analise == 0:
            return 'Тематическое моделирование'
        if self.type_analise == 1:
            return 'Факты EBM'
        if self.type_analise == 2:
            return 'Отрисовка графа'

    def get_status(self):
        return task_statuses[self.status][1]

    def time_delta(self):
        if self.end_date is None:
            return None
        delta_time = (self.end_date - self.start_date).seconds
        return delta_time