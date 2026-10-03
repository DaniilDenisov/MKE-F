FROM nginx:1.28-alpine

COPY docker/nginx.conf /etc/nginx/nginx.conf
COPY preprocessor /usr/share/nginx/html/preprocessor
COPY postprocessor /usr/share/nginx/html/postprocessor
COPY reference /usr/share/nginx/html/reference
COPY shared /usr/share/nginx/html/shared
COPY src /usr/share/nginx/html/src
COPY scripts /usr/share/nginx/html/scripts
COPY examples /usr/share/nginx/html/examples

ARG MKEF_VERSION=dev
ARG VCS_REF=unknown
LABEL org.opencontainers.image.title="MKE-F web" \
      org.opencontainers.image.description="Web interface and reference for MKE-F" \
      org.opencontainers.image.source="https://github.com/DaniilDenisov/MKE-F" \
      org.opencontainers.image.version="${MKEF_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

USER nginx
EXPOSE 8080
