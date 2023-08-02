from rest_framework import serializers
from rest_framework.authtoken.models import Token
from django.contrib.auth import get_user_model
from .models import UserPermissions
User = get_user_model()


class UserLoginSerilizer(serializers.Serializer):
    model = User

    email = serializers.EmailField(required=True)
    password = serializers.CharField(required=True)


class TokenSeriazliser(serializers.ModelSerializer):

    class Meta:
        model = Token
        fields = ['key']


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