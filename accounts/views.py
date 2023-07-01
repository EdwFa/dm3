from django.contrib.auth import authenticate, login, logout

from django.core.exceptions import ObjectDoesNotExist
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.request import Request
from rest_framework.permissions import AllowAny

from .serializers import *


class LoginApi(APIView):
    permission_classes = [AllowAny]

    def post(self, request: Request):
        """Настроить защиту позже"""
        print('Someone try to login...')
        serializer = UserLoginSerilizer(data=request.data)
        print(serializer)
        if serializer.is_valid():
            print(serializer.validated_data)
            authenticated_user = authenticate(**serializer.validated_data)
            if authenticated_user:
                print(f'Welcome {authenticated_user}')
                try:
                    token = Token.objects.get(user=authenticated_user)
                    token.delete()
                except ObjectDoesNotExist:
                    pass
                token = Token.objects.create(user=authenticated_user)
                print(token)
                return Response(TokenSeriazliser(token).data)
            return Response(serializer.errors, status=403)
        else:
            return Response(serializer.errors, status=400)


class LogoutApi(APIView):

    def get(self, request: Request):
        if request.user.is_authenticated and request.user.is_active:
            print(f'Logout {request.user.username}...')
            logout(request)
            return Response(data={'Status': 'Success'}, status=200)
        else:
            return Response({'Status': 'Unknowed user'}, status=403)




