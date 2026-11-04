FROM python:3.11.9-bookworm AS numeric
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
WORKDIR /app
COPY section_solver/requirements.txt /tmp/section-requirements.txt
RUN pip install --no-cache-dir -r /tmp/section-requirements.txt

FROM numeric AS runtime
COPY section_service/requirements.txt /tmp/service-requirements.txt
RUN pip install --no-cache-dir -r /tmp/service-requirements.txt
ENV HOME=/tmp MPLCONFIGDIR=/tmp/matplotlib
COPY section_solver /app/section_solver
COPY section_service /app/section_service
USER 10001:10001
EXPOSE 8000
CMD ["uvicorn","section_service.app:app","--host","0.0.0.0","--port","8000","--workers","1"]
