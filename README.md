# МКЭ-Ф

Учебная программа для расчёта плоских ферм и рам методом конечных элементов.

## Запуск в Docker

Установите Docker Desktop или Docker Engine с Compose.

### Готовый релиз из Docker Hub

После клонирования этого репозитория выполните в его корне:

```sh
docker compose -f compose.release.yaml up -d
```

Docker Compose автоматически загрузит и запустит образы релиза `0.7.1`:

- [denisovds/mkef-solver](https://hub.docker.com/r/denisovds/mkef-solver)
- [denisovds/mkef-web](https://hub.docker.com/r/denisovds/mkef-web)

Какой именно Git-коммит записан в образ, можно проверить по стандартной
OCI-метке:

```sh
docker image inspect denisovds/mkef-web:0.7.1 --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
docker image inspect denisovds/mkef-solver:0.7.1 --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

Новые релизные образы собираются скриптом `scripts/release_images.py`. Он
отказывается работать с незакоммиченными изменениями, сверяет версию с
`CHANGELOG.md` и `compose.release.yaml`, добавляет полный Git SHA в OCI-метку и
создаёт дополнительный идентифицирующий тег вида `git-<короткий SHA>`. Для
строгой фиксации конкретного содержимого образа следует использовать digest
`sha256:...`, который Docker показывает после публикации.

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
