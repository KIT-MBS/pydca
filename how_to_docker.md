# pyDCA: Docker Workaround

Create a Dockerfile with the following content:
```bash
FROM --platform=linux/amd64 python:3.7.8

RUN pip install pydca
WORKDIR /data
```

Build your docker image in the terminal:
```bash
$ docker build -t pydca:latest
```


Run a new container interactively in your terminal:
```bash
$ mkdir -p data
$ docker run -it --platform linux/amd64 --volume $(pwd)/data:/data pydca:latest bash
```
