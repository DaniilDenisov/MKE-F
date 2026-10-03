FROM nginx:1.28-alpine

COPY docker/nginx.conf /etc/nginx/nginx.conf
COPY preprocessor /usr/share/nginx/html/preprocessor
COPY postprocessor /usr/share/nginx/html/postprocessor
COPY reference /usr/share/nginx/html/reference
COPY shared /usr/share/nginx/html/shared
COPY src /usr/share/nginx/html/src
COPY scripts /usr/share/nginx/html/scripts
COPY examples /usr/share/nginx/html/examples

USER nginx
EXPOSE 8080
