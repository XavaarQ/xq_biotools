# syntax=docker/dockerfile:1
#artifactory.releng.ginkgo.zone/docker/
FROM artifactory.releng.ginkgo.zone/docker/python:3.11

# set a directory for the app
WORKDIR /usr/src/app

# copy all the files to the container
COPY . .

# install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# define the port number the container should expose
# EXPOSE 5000

# run the command
ENTRYPOINT ["python", "./codon_harmonizer.py"]
