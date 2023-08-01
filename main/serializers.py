from rest_framework import serializers
from .models import Article, TaskSearch, TaskAnalise
from accounts.models import UserPermissions


class ArticleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Article
        fields = ('uid', 'aid', 'titl', 'mesh', 'majr', 'subh', 'auth', 'jour', 'affl', 'pdat', 'tiab', 'ptyp', 'url', 'urlaid', 'pt', 'pl')


class TaskSearchSerializer(serializers.ModelSerializer):
    class Meta:
        model = TaskSearch
        fields = '__all__'


class TaskAnaliseSerializer(serializers.ModelSerializer):
    class Meta:
        model = TaskAnalise
        fields = '__all__'


class UserPermissionsSerialiser(serializers.ModelSerializer):
    class Meta:
        model = UserPermissions
        fields = '__all__'