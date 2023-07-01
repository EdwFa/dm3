from django.urls import path
from .views import *


urlpatterns = [
    path('', CreateTaskView.as_view(), name='index'),
    path('articles', CheckStatusTaskView.as_view(), name='check'),
    path('graphs', GetGraphData.as_view(), name='graph'),
]