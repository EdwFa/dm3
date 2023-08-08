FROM python:3.10.12

WORKDIR /usr/src/summariseapi

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN pip install --upgrade pip
RUN pip install pyproject-toml
RUN pip install python-dev-tools --user --upgrade
# RUN pip install flask
COPY ./requirements.txt requirements.txt
RUN pip install -r requirements.txt

COPY . .

EXPOSE 5000

CMD ["python", "main.py"]