# МКЭ-Ф

[for English - see below](#english)

Учебная программа для расчёта плоских ферм и рам методом конечных элементов.

## Работа с лабораторными

Docker-развёртывание удобно для расчётов через веб-интерфейс, но плохо подходит
для лабораторных: Octave работает внутри контейнера, а команды из справочника
предполагают интерактивную работу с исходным кодом и матрицами.

Для лабораторных установите GNU Octave непосредственно в свою систему,
запустите его и выберите папку репозитория `MKE-F` в качестве текущей рабочей
папки (в ней находится файл `setup.m`). Команды из лабораторных выполняйте
в командном окне Octave. Например, для запуска всех лабораторных:

```octave
setup;
addpath(fullfile(pwd, 'reference', 'examples'));
results = run_reference_examples();
```

Задания и пояснения находятся в [учебном справочнике](reference/RU/index.html).

## Запуск в Docker

Установите Docker Desktop или Docker Engine с Compose.

### Готовый релиз из Docker Hub

После клонирования этого репозитория выполните в его корне:

```sh
docker compose -f compose.release.yaml up -d
```

Docker Compose автоматически загрузит и запустит образы релиза `0.8.3`:

- [denisovds/mkef-solver](https://hub.docker.com/r/denisovds/mkef-solver)
- [denisovds/mkef-web](https://hub.docker.com/r/denisovds/mkef-web)

Какой именно Git-коммит записан в образ, можно проверить по стандартной
OCI-метке:

```sh
docker image inspect denisovds/mkef-web:0.8.3 --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
docker image inspect denisovds/mkef-solver:0.8.3 --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
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

Подробности находятся в [справочнике](reference/index.html):
[Русский](reference/RU/index.html) · [English](reference/EN/index.html).
Запуск, настройка и проверка проекта:
[Русский](reference/RU/project-operations.html) · [English](reference/EN/project-operations.html).

Порядок выпуска новой версии описан в [руководстве разработчика](DEVGUIDE.md).


### Узловые ограничения и MPC

Для рамы доступны все семь масок закреплений. Типы `5`, `6`, `7` закрепляют
только `ux`, `uy`, `thetaZ`. Ферма использует три уникальные комбинации.
Редактор сохраняет маски при смене семейства и требует убрать `thetaZ`
перед переходом к ферме.

Однородные MPC доступны в статике, модальном и переходном расчёте.
Редактор MPC и помощник «узел на оси» создают секцию `mpc`; помощник сохраняет
коэффициенты после правки геометрии. См. [формат и знаки сил](postprocessor/schema-v3.md).
Примеры: `CaseMPCFrame`, `CaseMPCModal`, `CaseMPCTransient`, `CaseMPCUniform`,
`CaseMPCTruss` в `examples/cases`.

Для сквозной проверки сначала выполните в Octave `setup; export_mpc_examples`,
затем `node scripts/test_mpc_browser.cjs` с установленным Playwright
(или заданным `MKEF_PLAYWRIGHT_MODULE`). `node scripts/generate-support-catalog.cjs --check`
проверяет синхронизацию общего справочника опор.

### Концевые освобождения Mz

Для каждого конца рамы 113 доступны независимые освобождения Mz во всех трёх анализах. Секция `releases` сохраняет прежний формат элементов. Внутренние повороты учитываются с исходной согласованной массой и доступны в карточках и временных графиках. Закрепление отсутствующего узлового поворота игнорируется с предупреждением; нагрузка Mz и ссылки MPC на него отклоняются.

См. [JSON v4](postprocessor/schema-v4.md) и примеры `CaseReleaseStatic`, `CaseReleaseModal`, `CaseReleaseTransient`. Проверка: `setup; addpath(fullfile(pwd, 'tests')); test_releases`, затем `node scripts/test_releases_browser.cjs` с Playwright.

### Линейная распределённая нагрузка

В статике для рам 113 доступны Uniform и Linear. В Linear задаются `qx1, qy1, qx2, qy2`
на концах элемента в местных или глобальных осях. Нулевой конец даёт треугольник,
разные значения — трапецию; допускается смена знака. Нагрузка действует на всю
длину элемента и суммируется с другими распределёнными и узловыми нагрузками.

Входной блок: `eload_linear`, число записей, строки
`21,elementId,coordinateSystem,qx1,qy1,qx2,qy2`.
См. [JSON v5](postprocessor/schema-v5.md), [формат](reference/RU/10-input-format.html#linear-loads)
и примеры `CaseLinearFrame`, `CaseTriangleFrame`, `CaseLinearRelease`, `CaseLinearMPC`.
Старый `eload_uniform` совместим без изменений. Эпюры N/V/M учитывают точное
линейное распределение и внутренние экстремумы.

Проверка: `setup; addpath('tests'); test_linear_loads` в Octave, затем
`node scripts/test_linear_loads_browser.cjs` с Playwright.

### NAFEMS Challenge Problem 5

Сохранённые входы: `examples/cases/nafems-challenge-5/manifest.csv`.
Короткая инструкция: [RU](reference/RU/06a-nafems-challenge-5.html) ·
[EN](reference/EN/06a-nafems-challenge-5.html).
Из корня проекта в Octave: `setup; run_nafems_challenge5('quick');`
или `run_nafems_challenge5('full');`. Результаты — в
`output/nafems-challenge-5/`. После quick браузерная проверка:
`node scripts/test_nafems_browser.cjs` (Playwright).
Исследование использует текущие рамные элементы и матрицы массы; расхождение
с опубликованными 22 формами указано в инструкции, параметры не подгоняются.

---

# English

MKE-F is an educational program for analyzing planar trusses and frames using the finite element method.

## Working with the labs

The Docker deployment is convenient for calculations through the web interface,
but is not well suited to the labs: Octave runs inside a container, while the
handbook commands assume interactive work with source code and matrices.

For the labs, install GNU Octave directly on your system, launch it, and select
the `MKE-F` repository folder as the current working directory (the folder
containing `setup.m`). Run the lab commands in the Octave command window.
For example, to run all labs:

```octave
setup;
addpath(fullfile(pwd, 'reference', 'examples'));
results = run_reference_examples();
```

Assignments and explanations are available in the [educational handbook](reference/EN/index.html).

## Running with Docker

Install Docker Desktop or Docker Engine with Compose.

### Prebuilt release from Docker Hub

After cloning this repository, run the following from its root directory:

```sh
docker compose -f compose.release.yaml up -d
```

Docker Compose will automatically download and start the `0.8.3` release images:

- [denisovds/mkef-solver](https://hub.docker.com/r/denisovds/mkef-solver)
- [denisovds/mkef-web](https://hub.docker.com/r/denisovds/mkef-web)

To check which Git commit an image contains, inspect its standard OCI label:

```sh
docker image inspect denisovds/mkef-web:0.8.3 --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
docker image inspect denisovds/mkef-solver:0.8.3 --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

New release images are built with `scripts/release_images.py`. The script refuses
to run with uncommitted changes, checks the version against `CHANGELOG.md` and
`compose.release.yaml`, adds the full Git SHA to an OCI label, and creates an
additional identifying tag in the form `git-<short SHA>`. To pin the exact image
contents, use the `sha256:...` digest reported by Docker after publishing.

Once the application starts, open <http://127.0.0.1:8080>.

To stop the application:

```sh
docker compose -f compose.release.yaml down
```

### Building from source

To build the images locally, run the following from the repository root:

```sh
docker compose up --build
```

If you use the standalone Compose command: `docker-compose up --build`.

Once the application starts, open <http://127.0.0.1:8080>.

Further details are available in the [handbook](reference/index.html):
[Русский](reference/RU/index.html) · [English](reference/EN/index.html).
For setup, configuration, and verification:
[Русский](reference/RU/project-operations.html) · [English](reference/EN/project-operations.html).

The release procedure is described in the [developer guide](DEVGUIDE.md) (in Russian).

### Nodal constraints and MPCs

All seven restraint masks are available for frames. Types `5`, `6`, and `7`
restrain only `ux`, `uy`, and `thetaZ`, respectively. Trusses use three unique
combinations. The editor preserves the masks when switching element families
and requires removing the `thetaZ` restraint before switching to trusses.

Homogeneous multipoint constraints (MPCs) are supported in static, modal, and
transient analyses. The MPC editor and the “node on axis” helper create the
`mpc` section; the helper preserves the coefficients after geometry edits.
See the [format and force sign conventions](postprocessor/schema-v3.md).
Examples: `CaseMPCFrame`, `CaseMPCModal`, `CaseMPCTransient`, `CaseMPCUniform`,
and `CaseMPCTruss` in `examples/cases`.

For an end-to-end check, first run `setup; export_mpc_examples` in Octave,
then run `node scripts/test_mpc_browser.cjs` with Playwright installed
(or `MKEF_PLAYWRIGHT_MODULE` set). Run `node scripts/generate-support-catalog.cjs --check`
to check that the shared support catalog is synchronized.

### Mz end releases

Independent Mz releases are available at each end of a type 113 frame element
in all three analysis types. The `releases` section preserves the existing
element format. Internal rotations are included using the original consistent
mass matrix and are available in result cards and time-history plots.
A restraint on an absent nodal rotation is ignored with a warning; an Mz load
or an MPC reference to that rotation is rejected.

See [JSON v4](postprocessor/schema-v4.md) and the `CaseReleaseStatic`,
`CaseReleaseModal`, and `CaseReleaseTransient` examples. To verify, run
`setup; addpath(fullfile(pwd, 'tests')); test_releases`, followed by
`node scripts/test_releases_browser.cjs` with Playwright.

### Linearly distributed loads

Static analysis of type 113 frame elements supports Uniform and Linear loads.
Linear loads specify `qx1, qy1, qx2, qy2` at the element ends in local or global
coordinates. A zero value at one end produces a triangular distribution;
different end values produce a trapezoidal distribution. Sign changes are
allowed. The load acts along the full element length and is added to other
distributed and nodal loads.

The input block consists of `eload_linear`, the number of records, and rows of
`21,elementId,coordinateSystem,qx1,qy1,qx2,qy2`.
See [JSON v5](postprocessor/schema-v5.md), the [input format](reference/EN/10-input-format.html#linear-loads),
and the `CaseLinearFrame`, `CaseTriangleFrame`, `CaseLinearRelease`, and
`CaseLinearMPC` examples. The existing `eload_uniform` format remains compatible
without changes. N/V/M diagrams account for the exact linear load distribution
and interior extrema.

To verify, run `setup; addpath('tests'); test_linear_loads` in Octave, followed by
`node scripts/test_linear_loads_browser.cjs` with Playwright.

### NAFEMS Challenge Problem 5

Saved input cases: `examples/cases/nafems-challenge-5/manifest.csv`.
Quick guide: [RU](reference/RU/06a-nafems-challenge-5.html) ·
[EN](reference/EN/06a-nafems-challenge-5.html).
From the project root in Octave, run `setup; run_nafems_challenge5('quick');`
or `run_nafems_challenge5('full');`. Results are written to
`output/nafems-challenge-5/`. After the quick run, check the browser views with
`node scripts/test_nafems_browser.cjs` (Playwright).
The study uses the current frame elements and mass matrices. The discrepancy
with the 22 published modes is documented in the guide; parameters are not
tuned to match them.
