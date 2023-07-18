from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.core.exceptions import ObjectDoesNotExist

from celery.result import AsyncResult
from django_celery_results.models import TaskResult

from dm.settings import RETMAX, MAX_COUNT
from .tasks import *



class BaseTaskView(APIView):
    taskModel = TaskSearch  # Начальная модель для поиска задачи
    files = ['search_ncbi', 'tematic_analise', 'clust_graph', 'heapmap', 'heirarchy', 'embeddings', 'info_graph']  # Все возмодные файлы для записи данных
    worker_func = parse_records
    label = 'data'
    retmax = 10000 # RETMAX

    def check_working_task(self, request, **kwargs):
        # Модуль проверки наличии уже запущенных запросов в поиске данному пользователю
        worked_tasks = self.taskModel.objects.filter(status=0, user=request.user, **kwargs).order_by('-start_date')
        if worked_tasks.count() != 0:  # Если пользователь имеет уже запущенные запросы выводим их ему
            try:
                return worked_tasks[0], TaskResult.objects.get(task_id=worked_tasks[0].task_id)
            except ObjectDoesNotExist:
                return worked_tasks[0], None
        return None, None

    def get_path_to_file(self, pk, file_name):
        path_to_file = os.path.join('datasets', str(pk), file_name)
        if not os.path.exists(path_to_file):
            raise Exception(f'Not found {path_to_file}')

        return path_to_file

    def create_task(self, **kwargs):
        # Создаем новую задачу
        task = self.taskModel.objects.create(**kwargs)
        return task

    def update_task(self, current_task, **kwargs):
        # обновляем состояние нашей задачи
        current_tasks = self.taskModel.objects.filter(id=current_task.id)
        current_tasks.update(**kwargs)

    def save_data(self, current_worker: TaskResult, pk):
        # Если наш запрос успешно обработался то сохраняем данные а если нет то файл пуст
        """
            Если наш запрос успешно обработался то сохраняем данные а если нет то файл пуст
            Формат data = {
                key1: value1, (List of dicts or None)
                ...
                keyN: valueN
            }
        """
        data = AsyncResult(current_worker.task_id, app=self.worker_func).get()
        for file_name in self.files:
            self.write_data(data[file_name], pk, file_name)
        return

    def write_data(self, data, pk, file_name):
        # Сохраняем наши изменения при условии что запрос прошел успешно
        f = open(self.get_path_to_file(pk, f'{file_name}.json'), 'w')
        json.dump(data, f)
        f.close()

    def create_data_response(self, pk):
        data = dict()
        for file_name in self.files:
            f = open(self.get_path_to_file(pk, f'{file_name}.json'), 'r')
            data[file_name] = json.load(f)
            f.close()
        return data

    def response_data(self, error_status, **kwargs):
        # Вывод наших данных, если запрос неудачный то вывводим ошибку б это
        return Response(data=kwargs, status=error_status)

    def get(self, request, **kwargs):
        current_task, current_worker = self.check_working_task(request, **kwargs)  # Получаем наш первый запущенный воркер по возрастанию даты запроса

        if current_task is None: # Если воркер отсутвует значит у пользователя сейчас свободна очередь запросов
            data = self.create_data_response(request.user.id)
            return self.response_data(200, data=data, message='Запросов нет')  # Выводим его последний запрос из базы

        if current_worker is None:  # Пользователь отправил запрос но обработчик не принял его
            current_task.delete()
            return self.response_data(500, message='Запрос завершен с ошибкой!', label=self.label)  # Выводим ошибку об отсутвии обработки запросов

        if current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED':
             return self.response_data(202, message=current_task.message, label=self.label)

        if current_worker.status == 'FAILURE':
            self.update_task(current_task, status=2, end_date=datetime.now(), message='Запрос завершен с ошибкой!')
            return self.response_data(500, message=current_task.message, label=self.label)

        if current_worker.status == 'SUCCESS':
            self.update_task(current_task, status=1, end_date=datetime.now(), message='Запрос успешно завершен!')
            self.save_data(current_worker, request.user.id)

        data = self.create_data_response(request.user.id)
        return self.response_data(200, data=data)


class SearchTaskView(BaseTaskView):
    taskModel = TaskSearch  # Начальная модель для поиска задачи
    files = ['search_ncbi']  #
    worker_func = parse_records
    label = 'search_ncbi'

    def get(self, request):
        current_task, current_worker = self.check_working_task(request)  # Получаем наш первый запущенный воркер по возрастанию даты запроса

        if current_task is None:  # Если воркер отсутвует значит у пользователя сейчас свободна очередь запросов
            data = self.create_data_response(request.user.id)
            task = TaskSearch.objects.filter(status=1, user=request.user).order_by('-end_date')
            if task.count() > 0:
                task = TaskSearchSerializer(task[0], many=False).data
            else:
                task = {'task': {
                    'full_query': '',
                    'translation_stack': '',
                    'query': '',
                    'count': 0,
                }}
            return self.response_data(200, data=data, message='Запросов нет', task=task)  # Выводим его последний запрос из базы

        if current_worker is None:  # Пользователь отправил запрос но обработчик не принял его
            task = TaskSearchSerializer(current_task, many=False).data
            current_task.delete()
            return self.response_data(500,
                                      message='Запрос завершен с ошибкой!',
                                      label=self.label,
                                      task=task)  # Выводим ошибку об отсутвии обработки запросов

        if current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED':
            return self.response_data(202, message=current_task.message, label=self.label, task=TaskSearchSerializer(current_task, many=False).data)

        if current_worker.status == 'FAILURE':
            self.update_task(current_task, status=2, end_date=datetime.now(), message='Запрос завершен с ошибкой!')
            return self.response_data(500, message=current_task.message, label=self.label, task=TaskSearchSerializer(current_task, many=False).data)

        if current_worker.status == 'SUCCESS':
            self.update_task(current_task, status=1, end_date=datetime.now(), message='Запрос успешно завершен!')
            self.save_data(current_worker, request.user.id)

        data = self.create_data_response(request.user.id)
        return self.response_data(200, data=data, task=TaskSearchSerializer(current_task, many=False).data)

    def post(self, request):
        current_task, current_worker = self.check_working_task(request)  # Получаем наш первый запущенный воркер по возрастанию даты запроса
        if (current_task is not None) and (current_worker is not None) and (current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED'):
            return self.response_data(403, message='В настоящее время вы не можете создать еще один запрос, дождитесь оканчания предыдущего.')

        query = create_query(**request.data)
        full_query, translation_stack, count = get_records(query)
        new_task = self.create_task(query=query, count=count, full_query=full_query, user=request.user,
                               translation_stack=translation_stack)
        task = parse_records.delay(query=query, count=count, new_task_id=new_task.id, retmax=self.retmax)
        new_task.message = 'Запрос получен'
        new_task.task_id = task.id
        new_task.save()
        data = {
            'data': TaskAnaliseSerializer(new_task, many=False).data,
        }
        return self.response_data(201, data=data)


class TematicAnaliseView(BaseTaskView):
    taskModel = TaskAnalise  # Начальная модель для поиска задачи
    files = ['tematic_analise', 'clust_graph', 'heapmap', 'heirarchy']
    worker_func = analise_records
    label = 'tematic_analise'

    def save_data(self, current_worker, id):
        return None

    def get(self, request):
        current_task, current_worker = self.check_working_task(request, type_analise=0)  # Получаем наш первый запущенный воркер по возрастанию даты запроса

        if current_task is None: # Если воркер отсутвует значит у пользователя сейчас свободна очередь запросов
            data = self.create_data_response(request.user.id)
            return self.response_data(200, data=data, message='Запросов нет')  # Выводим его последний запрос из базы

        if current_worker is None:  # Пользователь отправил запрос но обработчик не принял его
            current_task.delete()
            return self.response_data(500, message='Запрос завершен с ошибкой!', data={'tematic_review': None})  # Выводим ошибку об отсутвии обработки запросов

        if current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED':
             return self.response_data(202, message=current_task.message, data={'tematic_analise': None})

        if current_worker.status == 'FAILURE':
            self.update_task(current_task, status=2, end_date=datetime.now(), message='Запрос завершен с ошибкой!')
            return self.response_data(500, message=current_task.message)

        if current_worker.status == 'SUCCESS':
            self.update_task(current_task, status=1, end_date=datetime.now(), message='Запрос успешно завершен!')
            self.save_data(current_worker, request.user.id)

        data = self.create_data_response(request.user.id)
        return self.response_data(200, data=data)

    def post(self, request):
        current_task, current_worker = self.check_working_task(request)  # Получаем наш первый запущенный воркер по возрастанию даты запроса
        if (current_task is not None) and (current_worker is not None) and (current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED'):
            return self.response_data(403, message='В настоящее время вы не можете создать еще один запрос, дождитесь оканчания предыдущего.')

        if len(request.data['articles']) < 2:
            return self.response_data(400, message='В настоящее время вы не можете создать еще один запрос, дождитесь оканчания предыдущего.')

        new_task = self.create_task(user=request.user, type_analise=0)
        task = analise_records.delay(pk=request.user.id, IdList=request.data['articles'], new_task_id=new_task.id)

        new_task.task_id = task.id
        new_task.message = 'Начинаем обработку...'
        new_task.save()
        data = {
            'data': TaskAnaliseSerializer(new_task, many=False).data,
        }
        return self.response_data(200, data=data)


class EmbeddingTaskView(BaseTaskView):
    taskModel = TaskAnalise  # Начальная модель для поиска задачи
    files = ['embeddings']
    worker_func = get_ddi_articles
    label = 'data'

    def get(self, request):
        current_task, current_worker = self.check_working_task(request, type_analise=1)  # Получаем наш первый запущенный воркер по возрастанию даты запроса
        if current_task is None:
            return self.response_data(404, message='Сделайте запрос',
                                      label=self.label)  # Выводим ошибку об отсутвии обработки запросов
        if current_worker is None:  # Пользователь отправил запрос но обработчик не принял его
            current_task.delete()
            return self.response_data(500, message='Ваш запрос не обработался! Пожайлуста повторите попытку позже',
                                      label=self.label)  # Выводим ошибку об отсутвии обработки запросов

        if current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED':
            return self.response_data(202, message=current_task.message, data=None)

        if current_worker.status == 'FAILURE':
            self.update_task(current_task, status=2, end_date=datetime.now(), message='Запрос завершен с ошибкой!')
            return self.response_data(500, message=current_task.message, label=self.label)

        if current_worker.status == 'SUCCESS':
            self.update_task(current_task, status=1, end_date=datetime.now(), message='Запрос успешно завершен!')
            data = AsyncResult(current_worker.task_id, app=summarise_text).get()
        return self.response_data(200, data=data)

    def post(self, request):
        current_task, current_worker = self.check_working_task(request)  # Получаем наш первый запущенный воркер по возрастанию даты запроса
        if (current_task is not None) and (current_worker is not None) and (current_worker.status == 'PROGRESS' or current_worker.status == 'STARTED'):
            return self.response_data(403, message='В настоящее время вы не можете создать еще один запрос, дождитесь оканчания предыдущего.')

        new_task = self.create_task(user=request.user, type_analise=1)
        task = get_ddi_articles.delay(new_task_id=new_task.id, **request.data)

        new_task.task_id = task.id
        new_task.message = 'Начало обработки...'
        new_task.save()

        data = {
            'data': TaskAnaliseSerializer(new_task, many=False).data,
            'message': 'Начало обработки...'
        }
        return self.response_data(200, data=data)

    def delete(self, request):
        for file_name in self.files:
            f = open(get_path_to_file(request.user.id, f'{file_name}.json'), 'w')
            f.close()
        data = {
            'data': None,
            'message': "Успешно очищено"
        }
        return self.response_data(200, data=data)


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


class SummariseEmbApi(APIView):
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
                'data': AsyncResult(worker_id.task_id, app=summarise_emb).get()
            }
        return Response(data=data, status=status.HTTP_200_OK)

    def post(self, request):
        print(request.data)
        task = summarise_emb.delay(pk=request.user.id)
        data = {
            'data': task.id,
        }
        return Response(data=data, status=status.HTTP_200_OK)


class MartUpApi(APIView):
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
                'data': AsyncResult(worker_id.task_id, app=summarise_emb).get()
            }
        return Response(data=data, status=status.HTTP_200_OK)

    def post(self, request):
        print(request.data)
        task = markup_artcile.delay(record=request.data['article'])
        data = {
            'data': task.id,
        }
        return Response(data=data, status=status.HTTP_200_OK)

class GetGraphData(BaseTaskView):
    files = ['info_graph']

    def get(self, request):
        data = self.create_data_response(request.user.id)
        return self.response_data(200, data=data)

