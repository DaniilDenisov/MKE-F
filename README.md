# МКЭ-Ф

Учебная программа для расчёта плоских ферм и рам методом конечных элементов.

## Запуск в Docker

Установите Docker Desktop или Docker Engine с Compose.

### Готовый релиз из Docker Hub

После клонирования этого репозитория выполните в его корне:

```sh
docker compose -f compose.release.yaml up -d
```

Docker Compose автоматически загрузит и запустит образы релиза `0.7.0`:

- [denisovds/mkef-solver](https://hub.docker.com/r/denisovds/mkef-solver)
- [denisovds/mkef-web](https://hub.docker.com/r/denisovds/mkef-web)

После запуска откройте <http://127.0.0.1:8080>.

Остановить приложение:

```sh
docker compose -f compose.release.yaml down
```

### Сборка из исходников

Чтобы собрать образы локально, выполните в корне репозитория:

```sh
docker compose up --build
```

Если используется отдельная команда Compose: `docker-compose up --build`.

После запуска откройте <http://127.0.0.1:8080>.

Подробности находятся в [справочнике](reference/index.html), включая
[настройку и проверку проекта](reference/project-operations.html).
