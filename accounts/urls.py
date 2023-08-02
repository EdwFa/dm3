from django.urls import path
from .views import *


urlpatterns = [
    path('login', LoginApi.as_view(), name='login'),
    path('logout', LogoutApi.as_view(), name='logout'),
    path('create', CreateUserApi.as_view(), name='create'),
]