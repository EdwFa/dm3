from rest_framework import serializers
from .models import Article, Task


class ArticleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Article
        fields = ('uid', 'aid', 'titl', 'mesh', 'majr', 'subh', 'auth', 'jour', 'affl', 'pdat', 'tiab', 'ptyp', 'url', 'urlaid', 'pt', 'pl')


class TaskSerializer(serializers.ModelSerializer):
    class Meta:
        model = Task
        fields = '__all__'