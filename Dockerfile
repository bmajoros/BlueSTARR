ARG TENSORFLOW_VERSION=2.15.0

FROM tensorflow/tensorflow:${TENSORFLOW_VERSION}-gpu
ARG TENSORFLOW_VERSION

RUN apt-get update && apt-get install --no-install-recommends -y git

WORKDIR /bluestarr

COPY non-tensorflow-reqs.txt .

# We have to specify the tensorflow version or else keras_nlp will force an upgrade to the latest
# tensorflow version
RUN pip install --no-cache-dir -r non-tensorflow-reqs.txt tensorflow==${TENSORFLOW_VERSION}
COPY *.py .
