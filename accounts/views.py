from django.contrib.auth import authenticate, login, logout

from django.core.exceptions import ObjectDoesNotExist
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.request import Request
from rest_framework.permissions import AllowAny

import os
import json

from .serializers import *


needed_files = ['search_ncbi', 'tematic_analise', 'clust_graph', 'heapmap', 'heirarchy', 'embeddings', 'info_graph', 'info_graph_journals', 'info_graph_countries', 'DTM']


def get_path_to_file(pk, file_name):
    path_to_file = os.path.join('datasets', str(pk), file_name)
    if not os.path.exists(path_to_file):
        f = open(path_to_file, 'w')
        json.dump(None, f)
        f.close()


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

                if not os.path.exists(f'datasets/{authenticated_user.id}'):
                    os.mkdir(f'datasets/{authenticated_user.id}')
                for need_file in needed_files:
                    get_path_to_file(authenticated_user.id, f'{need_file}.json')

                print(token)
                return Response(data={'token': TokenSeriazliser(token).data,
                                      'email': authenticated_user.email,
                                      'allow': authenticated_user.allow_status}, status=200)
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




