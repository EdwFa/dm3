FROM python:3.10.12

WORKDIR /usr/src/datamed

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN pip install --upgrade pip
RUN pip install python-dev-tools --user --upgrade
COPY ./requirements.txt /usr/src/datamed/requirements.txt
RUN pip install -r requirements.txt

COPY . /usr/src/datamed/

EXPOSE 5000

CMD ['python', 'main.py']