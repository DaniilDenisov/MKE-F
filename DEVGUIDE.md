# Руководство разработчика

## Две языковые версии справочника

Страницы хранятся отдельно в `reference/RU/` и `reference/EN/`. При изменении
главы обновляйте обе версии, сохраняя одинаковые имена файлов и идентификаторы
разделов. Подписи SVG и строки интерактивных сценариев находятся в `assets/`
соответствующего языка. Общие стили и перенаправление старых ссылок находятся
в `reference/assets/`, исполняемые лабораторные — в `reference/examples/`.
Корневой `reference/index.html` выбирает язык; остальные старые HTML-адреса
перенаправляют в RU, сохраняя параметры и якорь ссылки.

Проверка структуры и лабораторных в GNU Octave из корня репозитория:

```octave
addpath(pwd);
setup;
addpath(fullfile(pwd, 'tests'));
test_reference_site;
test_reference_examples;
```

Проверка браузера требует установленного Playwright. При необходимости
`MKEF_PLAYWRIGHT_MODULE` задаёт путь к модулю. Команда
`node scripts/test_reference_browser.cjs` проверяет открытие с диска,
переключение языков, старые ссылки, совпадение расчётов и мобильную вёрстку.
Чтобы дополнительно проверить работающий Docker-стек, задайте
`MKEF_REFERENCE_BASE_URL=http://127.0.0.1:8080`. Снимки сохраняются в
`output/reference-browser/`. Публикация релиза для этих проверок не нужна.

## Выпуск новой версии

Релиз состоит из одного Git-коммита, аннотированного Git-тега и двух Docker-
образов. Полный SHA исходного коммита записывается в OCI-метку
`org.opencontainers.image.revision` каждого образа.

Не переиспользуйте уже опубликованный номер версии. Для исправления релиза
увеличивайте patch-версию: например, после `0.7.1` выпускайте `0.7.2`.

### 1. Подготовить номер версии

Замените старую версию на новую в следующих местах:

1. Добавьте новую верхнюю секцию `## [X.Y.Z] - YYYY-MM-DD` в `CHANGELOG.md`.
2. Обновите теги обоих образов в `compose.release.yaml`.
3. Обновите номер релиза и примеры команд `docker image inspect` в `README.md`.

Убедитесь, что старый номер остался только в истории changelog:

```powershell
rg -n "0\.7\.1|0\.7\.2" README.md CHANGELOG.md compose.release.yaml
```

### 2. Проверить проект до релиза

Соберите тестовый образ решателя и запустите интеграционные тесты:

```powershell
docker build --target test -t mkef-solver-test -f docker/solver.Dockerfile .
docker run --rm mkef-solver-test
```

Затем проверьте полный локальный стек:

```powershell
docker compose up --build -d
./scripts/compose_smoke_test.ps1
docker compose down
```

### 3. Закоммитить и отправить исходники

Перед сборкой релиза все изменения должны находиться в одном Git-коммите, а
рабочее дерево должно быть чистым:

```powershell
git status --short
git add -A
git diff --cached --check
git diff --cached
git commit -m "Prepare Docker release X.Y.Z"
git push origin master
git status
```

Перед `git commit` внимательно просмотрите вывод `git diff --cached`: в релиз
не должны попасть временные файлы, секреты или посторонние изменения.

Последний `git status` должен показать, что локальная ветка совпадает с
`origin/master`, а незакоммиченных файлов нет. Release-скрипт остановится, если
эти условия не выполнены.

### 4. Проверить план релиза

Запустите release-скрипт без сборки и публикации:

```powershell
python scripts/release_images.py --dry-run
```

Проверьте в выводе:

- номер версии;
- полный Git SHA;
- имена `denisovds/mkef-solver` и `denisovds/mkef-web`;
- теги `X.Y.Z`, `git-<короткий SHA>` и `latest`.

### 5. Собрать и опубликовать Docker-образы

Docker должен быть запущен, а пользователь должен быть авторизован в Docker
Hub под учётной записью с правом записи в репозитории `denisovds`:

```powershell
docker login
python scripts/release_images.py --push
```

Скрипт выполняет следующие проверки и действия:

1. Проверяет чистоту рабочего дерева.
2. Берёт версию из верхней секции `CHANGELOG.md`.
3. Проверяет эту версию в `compose.release.yaml`.
4. Проверяет, что `HEAD` уже отправлен в upstream.
5. Создаёт локальный аннотированный тег `vX.Y.Z`.
6. Собирает solver- и web-образы с OCI-метками версии и полного Git SHA.
7. Проверяет метки собранных образов.
8. Отправляет для каждого образа теги `X.Y.Z`, `git-<короткий SHA>` и
   `latest`.

Если загрузка оборвалась из-за сети, повторите ту же команду. Уже загруженные
слои будут использованы повторно.

### 6. Отправить Git-тег

После успешной публикации обоих Docker-образов отправьте созданный скриптом
Git-тег:

```powershell
git push origin vX.Y.Z
```

Не перемещайте уже опубликованный Git-тег на другой коммит. Если после релиза
обнаружена ошибка, исправьте её и выпустите следующую patch-версию.

### 7. Проверить опубликованный релиз

Получите образы из Docker Hub:

```powershell
docker pull denisovds/mkef-web:X.Y.Z
docker pull denisovds/mkef-solver:X.Y.Z
```

Сравните SHA в обоих образах с коммитом Git-тега:

```powershell
git rev-list -n 1 vX.Y.Z
docker image inspect denisovds/mkef-web:X.Y.Z --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
docker image inspect denisovds/mkef-solver:X.Y.Z --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

Все три команды должны вывести один и тот же полный Git SHA. Неизменяемый
идентификатор конкретного содержимого Docker-образа — digest `sha256:...`,
показываемый командой `docker pull`.

Дополнительно проверьте удалённый Git-тег:

```powershell
git ls-remote origin 'refs/tags/vX.Y.Z^{}'
```

### 8. Обновить запущенный релиз

После публикации загрузите новые образы и пересоздайте контейнеры:

```powershell
docker compose -f compose.release.yaml pull
docker compose -f compose.release.yaml up -d --force-recreate
docker compose -f compose.release.yaml ps
```

Откройте <http://127.0.0.1:8080> и выполните `Ctrl+F5`, чтобы браузер не
показывал старые статические файлы из кэша.

### Частые ошибки

`Working tree is not clean`

: Закоммитьте необходимые изменения либо уберите посторонние изменения из
  рабочего дерева. Не выпускайте образ из незакоммиченного состояния.

`HEAD ... is not the commit at origin/master`

: Сначала выполните `git push origin master`, затем повторите release-скрипт.

`Git tag vX.Y.Z points to ... not HEAD`

: Такой номер версии уже связан с другим коммитом. Увеличьте patch-версию; не
  перемещайте опубликованный тег.

`denied: requested access to the resource is denied`

: Выполните `docker login` под учётной записью с доступом к образам
  `denisovds/mkef-solver` и `denisovds/mkef-web`.

После обновления видны старые страницы

: Проверьте версию запущенных контейнеров, пересоздайте их командами из шага 8
  и обновите страницу через `Ctrl+F5`.
  Сообщение браузера `Unknown section marker "eload_uniform"` при новом
  решателе также может означать старую копию JavaScript в кэше.
  Nginx требует перепроверять статические файлы (`Cache-Control: no-cache`),
  но после перехода со старой версии нужна одна принудительная перезагрузка,
  чтобы браузер получил новую политику кэширования.
