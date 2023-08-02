from rest_framework import serializers
from .models import Article, TaskSearch, TaskAnalise, User
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
    topic = serializers.CharField(source='get_topic')
    class Meta:
        model = UserPermissions
        fields = '__all__'

class AdminPanelUsersSerialiser(serializers.ModelSerializer):
    permissions = UserPermissionsSerialiser(many=True)
    allow_status = serializers.CharField(source='get_allow')

    class Meta:
        model = User
        fields = ('email', 'is_admin', 'allow_status', 'permissions', 'last_login', 'task_search')

class AdminPanelSearchSerialiser(serializers.ModelSerializer):
    user = serializers.CharField(source='get_user')
    status = serializers.CharField(source='get_status')
    work_time = serializers.IntegerField(source='time_delta')

    class Meta:
        model = TaskSearch
        fields = ('user', 'query', 'full_query', 'translation_stack', 'start_date', 'end_date', 'status', 'count', 'work_time')

class AdminPanelAnaliseSerialiser(serializers.ModelSerializer):
    user = serializers.CharField(source='get_user')
    status = serializers.CharField(source='get_status')
    type_analise = serializers.CharField(source='get_type_analise')
    work_time = serializers.IntegerField(source='time_delta')

    class Meta:
        model = TaskAnalise
        fields = '__all__'
