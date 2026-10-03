FROM gnuoctave/octave:11.3.0 AS base

USER root

RUN apt-get update \
    && apt-get install --no-install-recommends --yes python3-venv \
    && rm -rf /var/lib/apt/lists/*

COPY backend/requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/mkef-venv \
    && /opt/mkef-venv/bin/pip install --no-cache-dir --requirement /tmp/requirements.txt \
    && rm /tmp/requirements.txt

RUN groupadd --gid 10001 mkef \
    && useradd --uid 10001 --gid 10001 --home-dir /tmp/mkef-home --no-create-home mkef \
    && mkdir -p /app /jobs /tmp/mkef-home \
    && chown -R mkef:mkef /jobs /tmp/mkef-home

WORKDIR /app
COPY --chown=mkef:mkef setup.m /app/setup.m
COPY --chown=mkef:mkef src /app/src
COPY --chown=mkef:mkef scripts /app/scripts
COPY --chown=mkef:mkef examples /app/examples
COPY --chown=mkef:mkef backend /app/backend

ENV PATH="/opt/mkef-venv/bin:${PATH}" \
    PYTHONPATH=/app \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    HOME=/tmp/mkef-home \
    MKEF_JOBS_ROOT=/jobs \
    MKEF_REPOSITORY_ROOT=/app \
    MKEF_RUNNER_SCRIPT=/app/scripts/run_solver_job.m

USER 10001:10001
EXPOSE 8000

FROM base AS test
USER root
COPY backend/requirements.txt backend/requirements-dev.txt /tmp/
RUN /opt/mkef-venv/bin/pip install --no-cache-dir --requirement /tmp/requirements-dev.txt \
    && rm /tmp/requirements.txt /tmp/requirements-dev.txt
COPY --chown=mkef:mkef tests /app/tests
COPY --chown=mkef:mkef README.md /app/README.md
COPY --chown=mkef:mkef .editorconfig /app/.editorconfig
COPY --chown=mkef:mkef reference /app/reference
COPY --chown=mkef:mkef preprocessor /app/preprocessor
COPY --chown=mkef:mkef postprocessor /app/postprocessor
COPY --chown=mkef:mkef shared /app/shared
USER 10001:10001
CMD ["pytest", "-q", "-p", "no:cacheprovider", "backend/tests"]

FROM base AS runtime
CMD ["uvicorn", "backend.app:app", "--host", "0.0.0.0", "--port", "8000", "--no-access-log"]
