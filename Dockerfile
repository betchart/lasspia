FROM python:2.7-slim

RUN mkdir /data
WORKDIR /lasspia
ADD . /lasspia
RUN pip install scipy astropy

VOLUME /data

CMD ["./lasspia.py","configs/cmassS.py","routines/quickscan.py"]
