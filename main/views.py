import json

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.core.exceptions import ObjectDoesNotExist

from django_celery_results.models import TaskResult
from celery.result import AsyncResult


from datetime import date, datetime
import time

from .tasks import *


class CreateTaskView(APIView):

    def get(self, request):

        retmax = 1000

        if check_working_task(request):
            return Response(data={'data': 'exist working tasks'}, status=status.HTTP_403_FORBIDDEN)

        query = create_query(request)
        full_query, translation_stack, count = get_records(query)
        new_task = create_task(query=query, count=count, full_query=full_query, user=request.user, translation_stack=translation_stack)

        print(query, count, new_task.id, full_query, translation_stack)
        task = parse_records.delay(query=query, count=count, new_task_id=new_task.id, retmax=retmax)

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
            with open('test_json.json') as f:
                data = json.load(f)
            # return Response(data={'data': data, 'message': 'No one query yet'}, status=status.HTTP_404_NOT_FOUND)
            return Response(data=data, status=status.HTTP_200_OK)

        try:
            worker_id = TaskResult.objects.get(task_id=current_task.task_id)
        except ObjectDoesNotExist:
            current_task.delete()
            return Response(data={'data': None, 'message': 'worker not found'}, status=status.HTTP_404_NOT_FOUND)

        if worker_id.status == 'PROGRESS' or worker_id.status == 'STARTED':
            return Response(data={'data': None, 'message': f'worker in progress. {current_task.message}'}, status=status.HTTP_202_ACCEPTED)

        if worker_id.status == 'FAILURE':
            current_task.status = 2
            current_task.end_date = datetime.now()
            current_task.save()
            return Response(data={'data': None, 'message': 'worker in failed'}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        if worker_id.status == 'SUCCESS':
            current_task.status = 1
            current_task.end_date = datetime.now()
            current_task.message = 'Done!'
            current_task.save()
            data = AsyncResult(worker_id.task_id, app=parse_records)
            with open('test_json.json', 'w') as f:
                json.dump(data.get(), f)

            return Response(data=data.get(), status=status.HTTP_200_OK)


class TematicAnaliseView(APIView):

    def get(self, request):
        try:
            current_task = TaskAnalise.objects.get(status=0, type_analise=0, user=request.user)
        except ObjectDoesNotExist:
            with open('test_analise_json.json') as f:
                table = json.load(f)
            with open('test_clust_graph.json') as f:
                graph = json.load(f)
            with open('test_heapmap.json') as f:
                heapmap = json.load(f)
            with open('test_heirarchy.json') as f:
                heirarchy = json.load(f)
            data = {
                'data': table,
                'graph': graph,
                'heapmap': heapmap,
                'heirarchy': heirarchy,
                'message': 'last query get'
            }
            return Response(data=data, status=status.HTTP_200_OK)

        print(current_task.task_id)
        try:
            worker_id = TaskResult.objects.get(task_id=current_task.task_id)
        except ObjectDoesNotExist:
            current_task.delete()
            return Response(data={'data': None, 'message': 'worker not found'}, status=status.HTTP_404_NOT_FOUND)

        if worker_id.status == 'PROGRESS' or worker_id.status == 'STARTED':
            return Response(data={'data': None, 'message': f'worker in progress. {current_task.message}'}, status=status.HTTP_202_ACCEPTED)

        if worker_id.status == 'FAILURE':
            current_task.status = 2
            current_task.end_date = datetime.now()
            current_task.save()
            return Response(data={'data': None, 'message': 'worker in failed'}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        if worker_id.status == 'SUCCESS':
            current_task.status = 1
            current_task.end_date = datetime.now()
            current_task.save()

        with open('test_analise_json.json') as f:
            table = json.load(f)
        with open('test_clust_graph.json') as f:
            graph = json.load(f)
        with open('test_heapmap.json') as f:
            heapmap = json.load(f)
        with open('test_heirarchy.json') as f:
            heirarchy = json.load(f)
        data = {
            'data': table,
            'graph': graph,
            'heapmap': heapmap,
            'heirarchy': heirarchy,
            'message': 'Done!'
        }
        return Response(data=data, status=status.HTTP_200_OK)

    def post(self, request):
        print(request.data)
        new_task = create_analise_task(user=request.user)
        if new_task is None:
            return Response(data={'data': None, 'message': 'worker in progress'}, status=status.HTTP_403_FORBIDDEN)
        task = analise_records.delay(IdList=request.data['articles'], new_task_id=new_task.id)

        new_task.task_id = task.id
        new_task.save()
        data = {
            'data': TaskAnaliseSerializer(new_task, many=False).data,
        }
        return Response(data=data, status=status.HTTP_200_OK)


class SummariseTextApi(APIView):
    def get(self, requset):
        task_id = requset.GET.get('task_id')
        try:
            worker_id = TaskResult.objects.get(task_id=task_id)
        except ObjectDoesNotExist:
            return Response(data={'data': None, 'message': 'worker not found'}, status=status.HTTP_404_NOT_FOUND)

        if worker_id.status == 'PROGRESS' or worker_id.status == 'STARTED':
            return Response(data={'data': None, 'message': 'worker in progress'}, status=status.HTTP_202_ACCEPTED)

        if worker_id.status == 'FAILURE':
            return Response(data={'data': None, 'message': 'worker in failed'}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        if worker_id.status == 'SUCCESS':
            data = {
                'data': AsyncResult(worker_id.task_id, app=summarise_text).get()
            }
        return Response(data=data, status=status.HTTP_200_OK)

    def post(self, request):
        print(request.data)
        task = summarise_text.delay(records=request.data['articles'])
        data = {
            'data': task.id,
        }
        return Response(data=data, status=status.HTTP_200_OK)

class DDIReviewApi(APIView):
    def get(self, request):

        try:
            current_task = TaskAnalise.objects.get(status=0, type_analise=1, user=request.user)
        except ObjectDoesNotExist:
            with open('test_ddi.json', 'r') as f:
                data = json.load(f)
            return Response(data=data, status=status.HTTP_200_OK)

        try:
            worker_id = TaskResult.objects.get(task_id=current_task.task_id)
        except ObjectDoesNotExist:
            current_task.delete()
            return Response(data={'data': None, 'message': 'worker not found'}, status=status.HTTP_404_NOT_FOUND)

        if worker_id.status == 'PROGRESS' or worker_id.status == 'STARTED':
            return Response(data={'data': None, 'message': f'worker in progress. {current_task.message}'}, status=status.HTTP_202_ACCEPTED)

        if worker_id.status == 'FAILURE':
            return Response(data={'data': None, 'message': 'worker in failed'}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        if worker_id.status == 'SUCCESS':
            current_task.status = 1
            current_task.end_date = datetime.now()
            current_task.save()
            with open('test_ddi.json', 'r') as f:
                try:
                    data = json.load(f)
                    data = data + AsyncResult(worker_id.task_id, app=get_ddi_articles).get()
                except:
                    data = AsyncResult(worker_id.task_id, app=get_ddi_articles).get()
            with open('test_ddi.json', 'w') as f:
                data = {
                    'data': data,
                    'message': 'Done!'
                }
                json.dump(data, f)

        with open('test_ddi.json', 'r') as f:
            data = json.load(f)
        return Response(data=data, status=status.HTTP_200_OK)

    def post(self, request):

        new_task = create_analise_task(user=request.user, type_analise=1)
        if new_task is None:
            return Response(data={'data': None, 'message': 'worker in progress'}, status=status.HTTP_403_FORBIDDEN)

        print(request.data)
        task = get_ddi_articles.delay(query=request.data['query'], new_task_id=new_task.id)

        new_task.task_id = task.id
        new_task.save()

        data = {
            'data': TaskAnaliseSerializer(new_task, many=False).data
        }
        return Response(data=data, status=status.HTTP_200_OK)

    def delete(self, request):
        f = open('test_ddi.json', 'w')
        f.close()
        return Response(data={'data': None, 'message': 'file is clear'}, status=status.HTTP_200_OK)

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
            nodes, edges = get_uniq_info_for_graph(articles[:50], 5, 0)
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

