from rest_framework import serializers
from django.contrib.auth.models import User
from rest_framework.authtoken.models import Token


class UserLoginSerilizer(serializers.Serializer):
    model = User

    username = serializers.CharField(required=True)
    password = serializers.CharField(required=True)

class TokenSeriazliser(serializers.ModelSerializer):

    class Meta:
        model = Token
        fields = ['key']