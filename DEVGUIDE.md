# Руководство разработчика

## Проверка операторов восстановления результатов

Экспортёр, HTTP runner, backend и просмотрщик работают только с v8 для 2D/3D.
Старые JSON v1–v7 нужно пересчитать из TXT. `formatVersion=8` можно указать явно;
опция `legacy` и любые другие версии отклоняются. Конвертера нет.
[Архитектура результатов](reference/result-recovery-architecture.md) и
[контракт v8](postprocessor/schema-v8.md).

В v8 режим `operatorStorage='auto'` объединяет строго одинаковые пары K/M
в `model.recoveryOperators`, а элементы хранят `recoveryOperatorId`.
Таблица используется только при уменьшении компактного JSON. Для явного выбора
доступны `inline` и `table`, например
`struct('formatVersion',8,'operatorStorage','inline')`. Чтение обоих вариантов
поддерживается; статические результаты матриц не содержат.

Все проверки выполнять в Docker. Единая команда из корня проекта:

```powershell
./scripts/test_recovery_pipeline.ps1
```

Нужны образы `mkef-solver-test:latest` и `mkef-browser-tests:latest`, собранные
из `docker/solver.Dockerfile` (target `test`) и `docker/browser-tests.Dockerfile`.
Образы предоставляют зависимости; текущие исходники монтируются только для
чтения. Каждый запуск получает новый каталог `output/recovery-runs/<run-id>`;
обычный `output` и эталоны предыдущих запусков не используются как вход.
При успехе появляется `recovery/PASS.json`, замеры — в
`recovery/browser-measurements.json`. Ошибка любого этапа останавливает запуск.
Замеры трёх вариантов хранения находятся в `recovery-pool/measurements.json`
и `recovery-pool/browser-measurements.json`: размер JSON/gzip, время экспорта,
чтения и память загруженных данных. Методика и решение по результатам —
в [отчёте](reference/recovery-operator-pooling-measurements.md).

Численные сценарии используют свежие v8-данные. Архивные JSON и отдельные
проверки v1–v7 удалены. Просмотрщик проверяет поддерживаемый формат и сообщает:
`Unsupported format. Recalculate.` Конвейер также запускает браузерные модульные
тесты и полный `backend/tests` с реальным Octave.
Новые примеры 2D, результаты и ссылки на решатель лежат в `planar-recovery/`
каталога запуска. `test_planar_recovery.m` и `test_planar_recovery.cjs` проверяют
K/M, кинематику, коэффициенты N/V/M, сечения, частичные истории и браузер.
Модуль `legacy-spatial-operators.js` удалён; локальные K/M берутся из решателя.
Редакторские формулы характеристик сечений и диагностики
ограничений намеренно сохраняются в обоих языках и требуют межъязыковых тестов.

`./scripts/test_v8_stack.ps1` собирает текущие runtime-образы и проверяет восемь
сценариев 2D/3D через nginx/API, Octave, браузер и повторное открытие скачанного
JSON и четыре готовых примера. Стек имеет уникальное имя,
не публикует порты и удаляется после проверки; рабочий сервер не затрагивается.
Готовые примеры пересчитываются в Docker функцией
`addpath('scripts'); export_postprocessor_examples();` из исходных TXT справочника.
Файлы из `output/postprocessor-examples` копируются в `examples/results`.

Назначенные сечения v8 содержат `elements[].sectionProperties` из модели
решателя. Геометрия каталога используется для контуров и сервиса Мизеса;
статический просмотрщик читает `sectionStress` вместе с критическими позициями.
Редактор продолжает вычислять характеристики для предварительного просмотра,
динамический постпроцессор — экстремумы текущего кадра. `stressField` вычисляет
поле в точке без поиска экстремумов; `MKEFSpatialSections.stress` сохраняет
поиск для динамики. Неиспользуемый `MKEFSections.stress` удалён из приложения;
его независимая формула находится только в `scripts/planar-stress-reference.cjs`.

`test_section_contract.m` экспортирует 42 случая в `section-contract/`.
`test_section_contract.cjs` сравнивает свойства/экстремумы с независимыми
формулами, проверяет 24 повреждённых контракта и отключает браузерные формулы
при загрузке результатов. `test_section_contract_browser.cjs` проверяет
окраску, сечения, SVG/CSV и тела запросов Мизеса с теми же запретами.
Оба набора включены в основной Docker-конвейер. Старые экспериментальные v8
без обязательных свойств/экстремумов нужно пересчитать из исходного TXT.

Границы общих модулей явные: `MKEFSections.geometry/contours` строят геометрию,
`properties` вычисляет характеристики для редактора, `resultProperties` читает
снимок решателя. У `MKEFSpatialSections` редактор использует `catalog`, а
просмотрщик — `resultCatalog`. Рисование сечения результата вызывает
`MKEFSectionDrawing.previewResult` с обязательными готовыми свойствами.
Не добавляйте выбор расчётного пути по наличию `model.recovery` или запасной
пересчёт характеристик при отсутствии полей v8.

`encodeResultJson.m` сохраняет конечные double через `%.17g`, нормализуя только
знак нуля. Нативный Octave `jsonencode` остаётся для строк и логических полей.
`test_result_numbers.m` / `.cjs` проверяют точную передачу 887 чисел, включая
очень малые K/M и истории. Эталоны берутся до сериализации, не из jsondecode.
`test_result_resources.m` / `.cjs` добавляют реальные расчёты 2D/3D и замеры
в `result-resources/`; оба набора входят в основной конвейер (14 JS-наборов).
Границы вычислений, обязательные парные проверки, замеры и пределы ресурсов
описаны в [архитектуре результатов](reference/result-recovery-architecture.md).

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

Релиз состоит из одного Git-коммита, аннотированного Git-тега и трёх Docker-
образов. Полный SHA исходного коммита записывается в OCI-метку
`org.opencontainers.image.revision` каждого образа.

Не переиспользуйте уже опубликованный номер версии. Для исправления релиза
увеличивайте patch-версию: например, после `0.7.1` выпускайте `0.7.2`.

### 1. Подготовить номер версии

Замените старую версию на новую в следующих местах:

1. Добавьте новую верхнюю секцию `## [X.Y.Z] - YYYY-MM-DD` в `CHANGELOG.md`.
2. Обновите теги всех трёх образов в `compose.release.yaml`.
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
- имена `denisovds/mkef-solver`, `denisovds/mkef-section-solver` и `denisovds/mkef-web`;
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
6. Собирает solver-, section-solver- и web-образы с OCI-метками версии и полного Git SHA.
7. Проверяет метки собранных образов.
8. Отправляет для каждого образа теги `X.Y.Z`, `git-<короткий SHA>` и
   `latest`.

Если загрузка оборвалась из-за сети, повторите ту же команду. Уже загруженные
слои будут использованы повторно.

### 6. Отправить Git-тег

После успешной публикации всех трёх Docker-образов отправьте созданный скриптом
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
docker pull denisovds/mkef-section-solver:X.Y.Z
```

Сравните SHA во всех трёх образах с коммитом Git-тега:

```powershell
git rev-list -n 1 vX.Y.Z
docker image inspect denisovds/mkef-web:X.Y.Z --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
docker image inspect denisovds/mkef-solver:X.Y.Z --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
docker image inspect denisovds/mkef-section-solver:X.Y.Z --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

Все четыре команды должны вывести один и тот же полный Git SHA. Неизменяемый
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
  `denisovds/mkef-solver`, `denisovds/mkef-section-solver` и `denisovds/mkef-web`.

После обновления видны старые страницы

: Проверьте версию запущенных контейнеров, пересоздайте их командами из шага 8
  и обновите страницу через `Ctrl+F5`.
  Сообщение браузера `Unknown section marker "eload_uniform"` при новом
  решателе также может означать старую копию JavaScript в кэше.
  Nginx требует перепроверять статические файлы (`Cache-Control: no-cache`),
  но после перехода со старой версии нужна одна принудительная перезагрузка,
  чтобы браузер получил новую политику кэширования.
