from django.urls import path
from .views import *


urlpatterns = [
    path('search/', SearchTaskView.as_view(), name='search'),
    path('analise/', TematicAnaliseView.as_view(), name='analise'),
    path('summarise', SummariseTextApi.as_view(), name='summarise'),

    path('graphs/', GetGraphData.as_view(), name='graph'),

    path('ddi_review', EmbeddingTaskView.as_view(), name='ddi'),
    path('summarise_emb', SummariseEmbApi.as_view(), name='summarise'),
    path('markup', MartUpApi.as_view(), name='markup'),

    path('download_vectors', DownloadVectors.as_view(), name='download_vectors'),
    path('download_metadata', DownloadMetaData.as_view(), name='download_metadata'),

    path('permissions', GetPermissions.as_view(), name='permissions'),
    path('translate', TranslateQuery.as_view(), name='translate'),

    path('admin', AdminPanelApi.as_view(), name='admin'),
]