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
ARG MKEF_VERSION=dev
ARG VCS_REF=unknown
LABEL org.opencontainers.image.title="MKE-F section solver" \
      org.opencontainers.image.description="Section FEM stress analysis API for MKE-F" \
      org.opencontainers.image.source="https://github.com/DaniilDenisov/MKE-F" \
      org.opencontainers.image.version="${MKEF_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"
USER 10001:10001
EXPOSE 8000
CMD ["uvicorn","section_service.app:app","--host","0.0.0.0","--port","8000","--workers","1"]
