from rest_framework import serializers
from rest_framework.authtoken.models import Token
from django.contrib.auth import get_user_model
User = get_user_model()


class UserLoginSerilizer(serializers.Serializer):
    model = User

    email = serializers.EmailField(required=True)
    password = serializers.CharField(required=True)

class TokenSeriazliser(serializers.ModelSerializer):

    class Meta:
        model = Token
        fields = ['key']