from django.urls import path
from .views import *


urlpatterns = [
    path('', CreateTaskView.as_view(), name='index'),
    path('articles/', CheckStatusTaskView.as_view(), name='check'),
    path('analise/', TematicAnaliseView.as_view(), name='analise'),
    path('graphs/', GetGraphData.as_view(), name='graph'),
    path('summarise', SummariseTextApi.as_view(), name='summarise'),

]