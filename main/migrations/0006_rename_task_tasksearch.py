# Generated by Django 4.2.2 on 2023-07-12 21:36

from django.conf import settings
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('main', '0005_task_message_taskanalise_message'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='Task',
            new_name='TaskSearch',
        ),
    ]
