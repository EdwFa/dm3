from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.request import Request
from rest_framework.decorators import api_view, permission_classes, authentication_classes
from rest_framework.permissions import IsAuthenticated, AllowAny
from rest_framework.authentication import TokenAuthentication, SessionAuthentication

from django_celery_results.models import TaskResult


from datetime import date, datetime
import time

from .serializers import *
from .tasks import *


class CreateTaskView(APIView):

    def get(self, request):
        if check_working_task(request):
            return Response(data={'data': 'exist working tasks'}, status=status.HTTP_403_FORBIDDEN)

        query = create_query(request)
        full_query, translation_stack, count = get_records(query)
        new_task = create_task(query=query, count=count, full_query=full_query, user=request.user, translation_stack=translation_stack)

        print(query, count, new_task.id, full_query, translation_stack)
        task = parse_records.delay(query=query, count=count, new_task_id=new_task.id)

        new_task.task_id = task.id
        new_task.save()

        data = {
            'task': TaskSerializer(new_task, many=False).data
        }

        return Response(data=data, status=status.HTTP_200_OK)

class CheckStatusTaskView(APIView):

    def get(self, request):
        # Получаем информацию о выполении задачи селери
        try:
            current_task = Task.objects.get(status=0, user=request.user)
        except ObjectDoesNotExist:
            last_task = Task.objects.filter(user=request.user).order_by('-end_date')
            if last_task.count() == 0:
                return Response(data={'data': None, 'message': 'No one query yet'}, status=status.HTTP_404_NOT_FOUND)

            last_task = last_task[0]
            start_time = time.time()
            articles = last_task.articles.all()
            print(articles.count())
            print("Time for get in db articles = ", time.time() - start_time)

            data = {
                'data': ArticleSerializer(articles, many=True).data,
                'columns': [
                    {'field': 'uid'},
                    {'field': 'titl'},
                    {'field': 'pdat'},
                    {'field': 'auth'},
                    {'field': 'jour'},
                    {'field': 'pt'}
                ],
                'message': 'no one task in progress'
            }
            start_time = time.time()
            print('Start get graph data...')
            # nodes, edges = get_uniq_info_for_graph(articles)
            # data['nodes'] = nodes
            # data['edges'] = edges
            print("Time for serialize articles = ", time.time() - start_time)
            return Response(data=data, status=status.HTTP_200_OK)
        try:
            worker_id = TaskResult.objects.get(task_id=current_task.task_id)
        except ObjectDoesNotExist:
            current_task.delete()
            return Response(data={'data': None, 'message': 'worker not found'}, status=status.HTTP_404_NOT_FOUND)

        if worker_id.status == 'PROGRESS' or worker_id.status == 'STARTED':
            return Response(data={'data': None, 'message': 'worker in progress'}, status=status.HTTP_403_FORBIDDEN)

        if worker_id.status == 'FAILURE':
            current_task.status = 2
            current_task.end_date = datetime.now()
            current_task.save()
            return Response(data={'data': None, 'message': 'worker in failed'}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        if worker_id.status == 'SUCCESS':
            current_task.status = 1
            current_task.end_date = datetime.now()
            current_task.save()
            data = {
                'data': ArticleSerializer(current_task.articles.all(), many=True).data,
                'columns': [
                    {'field': 'uid'},
                    {'field': 'titl'},
                    {'field': 'pdat'},
                    {'field': 'auth'},
                    {'field': 'jour'},
                    {'field': 'pt'}
                ],
                'message': 'worker is done success'
            }

            return Response(data=data, status=status.HTTP_200_OK)

class GetGraphData(APIView):

    def get(self, request):
        # Получаем данные для графа
        try:
            current_task = Task.objects.get(status=0, user=request.user)
            return Response(data={'data': None, 'message': 'worker in progress'}, status=status.HTTP_403_FORBIDDEN)
        except ObjectDoesNotExist:
            last_task = Task.objects.filter(user=request.user).order_by('-end_date')
            if last_task.count() == 0:
                return Response(data={'data': None, 'message': 'No one query yet'}, status=status.HTTP_404_NOT_FOUND)

            last_task = last_task[0]
            start_time = time.time()
            articles = last_task.articles.all()
            print(articles.count())
            print("Time for get in db articles = ", time.time() - start_time)

            start_time = time.time()
            print('Start get graph data...')
            nodes, edges = get_uniq_info_for_graph(articles[:100], 5)
            print("Time for create graph data = ", time.time() - start_time)
            start_time = time.time()
            data = {
                'graph': {
                    'nodes': nodes,
                    'edges': edges
                }
            }
            print("Time for serialize articles = ", time.time() - start_time)
            return Response(data=data, status=status.HTTP_200_OK)

