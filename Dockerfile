FROM python:3.10
ENV PYTHONUNBUFFERED 1
RUN python -m pip install --upgrade pip && pip --version
RUN mkdir /calc_vigas
WORKDIR /calc_vigas
COPY requirements.txt /calc_vigas/
COPY ./ /calc_vigas/
RUN pip install -r requirements.txt && ls -la && chmod -R 777 /calc_vigas/setup.sh
EXPOSE 80
ENTRYPOINT [ "/calc_vigas/setup.sh" ]
CMD bash
