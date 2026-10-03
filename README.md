# МКЭ-Ф

Учебная программа для расчёта плоских ферм и рам методом конечных элементов.

## Запуск через Docker

Установите Docker Desktop или Docker Engine с Compose и выполните в корне
репозитория:

```sh
docker compose up --build
```

Если используется отдельная команда Compose: `docker-compose up --build`.

После запуска откройте <http://127.0.0.1:8080>.

### Готовые образы Docker Hub

- [denisovds/mkef-solver](https://hub.docker.com/r/denisovds/mkef-solver)
- [denisovds/mkef-web](https://hub.docker.com/r/denisovds/mkef-web)

Загрузить образы релиза `0.7.0`:

```sh
docker pull denisovds/mkef-solver:0.7.0
docker pull denisovds/mkef-web:0.7.0
```

Для обоих образов также опубликован тег `latest`.

Подробности находятся в [справочнике](reference/index.html), включая
[настройку и проверку проекта](reference/project-operations.html).
