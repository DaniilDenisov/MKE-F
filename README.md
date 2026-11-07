# МКЭ-Ф

Учебная программа для расчёта плоских и пространственных ферм и рам методом конечных элементов.

[Сайт проекта](https://mkefweb.com/) · [English](#english)

## Возможности

- Статический, модальный и переходный расчёты на GNU Octave.
- Браузерный редактор 2D/3D: геометрия, нагрузки, опоры, многоточечные связи (MPC), концевые освобождения, каталог сечений и Undo/Redo.
- Единая WebGL-сцена: деформации, эпюры, реакции, собственные формы, анимация и графики. Нормальные напряжения и окраска статических результатов по максимальному |σx| от продольной силы и изгиба.
- Сохранение результатов в JSON v8, сцен 2D — в PNG/SVG, сцен 3D — в PNG; SVG графиков и сечений, CSV численных данных. Просмотр сохранённых JSON v8 без решателя; файлы v1–v7 требуют пересчёта исходного TXT; справочник RU/EN.

Входной TXT содержит обязательный блок `analysis`. Узловая нагрузка 10 — статическая; в transient используются гармоника 11, импульс 12 при t=dt и ступень 13.

Пространственные рамы принимают круг, прямоугольник и круглые, квадратные и прямоугольные трубы при импорте, редактировании, экспорте и расчёте. Депланация и смещение центра сдвига не учитываются.

В статическом постпроцессоре выберите стержень и координату: панель **Section stresses** покажет карты нормальных напряжений и Мизеса, с развёрнутым просмотром и экспортом двух карт в PNG. Сдвиг и кручение рассчитывает Python-сервис `sectionproperties` в Docker Compose; фермы и чистое растяжение/изгиб отображаются без сервиса. Панель **Member von Mises maxima · static** сохраняет максимумы по стержням, а [диагностическая страница](reference/section-fem-lab.html) — независимый ввод геометрии и усилий. [Просмотр и ограничения](reference/RU/08-result-recovery.html#section-view-title).

## Запуск

Установите Docker Desktop или Docker Engine с Compose. Для текущих возможностей выполните в корне клонированного репозитория:

```sh
docker compose up --build -d
```

Откройте [http://127.0.0.1:8080](http://127.0.0.1:8080). Остановка: `docker compose down`.

В препроцессоре список **Load example…** рядом со ссылкой **Examples** загружает готовую модель: по два простых примера для статического, модального и переходного расчётов в 2D и 3D (всего 12). Выберите пример и нажмите **Run analysis**. Примеры также доступны при открытии редактора без сервера; модель можно изменить и сохранить через **Download TXT**.

Для готового релиза:

```sh
docker compose -f compose.release.yaml up -d
```

[compose.release.yaml](compose.release.yaml) закреплён на версии **0.8.6** для всех трёх сервисов: решателя, сервиса сечений и веб-интерфейса. Остановка: `docker compose -f compose.release.yaml down`.

## Проверка перед релизом

В PowerShell 7+ при работающем Docker выполните
`pwsh -NoProfile -File ./scripts/testing/test_release.ps1`. Сценарий собирает образы,
запускает полную приёмку и сохраняет отчёт в
`output/release-runs/<run-id>/report.md`. [Порядок и параметры](DEVGUIDE.md#полная-приёмка-перед-релизом).

[Служебные команды](scripts/README.md) · [Структура и регистрация тестов](tests/README.md).

## Документация и лабораторные

- [Справочник](reference/RU/index.html) · [3D: быстрый старт](reference/RU/3d-quick-start.html) · [Пример модели](examples/cases/CaseSpatialQuickStart.txt).
- [Запуск, настройка и проверки](reference/RU/project-operations.html) · [История изменений](CHANGELOG.md).
- [Проверочный пример NAFEMS Challenge 5](reference/RU/06a-nafems-challenge-5.html).

Для лабораторных установите GNU Octave локально, откройте папку репозитория, выполните `setup` и следуйте заданиям справочника. Лабораторные предполагают интерактивную работу с кодом и матрицами; Docker предназначен прежде всего для веб-интерфейса.

---

# English

MKE-F is an educational finite element program for planar and spatial trusses and frames.

[Project website](https://mkefweb.com/)

## Features

- Static, modal and transient analysis using GNU Octave.
- Browser-based 2D/3D editor: geometry, loads, supports, multipoint constraints (MPC), end releases, section catalog and undo/redo.
- Shared WebGL scene: deformation, force diagrams, reactions, mode shapes, animation and plots. Normal stresses and static result coloring by maximum |σx| from axial force and bending.
- Save results as JSON v8, 2D scenes as PNG/SVG and 3D scenes as PNG; SVG charts and sections, CSV numerical data. View saved JSON v8 without the solver; v1–v7 files require recalculating the original TXT; RU/EN handbook.

Each input TXT requires an `analysis` block. Nodal load 10 is static only; transient analysis accepts harmonic load 11, pulse 12 at t=dt and step 13.

Spatial frames support solid circles, rectangles, circular tubes, square tubes and rectangular tubes at import, editing, export and calculation. Restrained warping and shear-center offsets are not modeled.

Select a member and coordinate in the static postprocessor: **Section stresses** shows normal and von Mises maps, with an expanded view and a combined PNG export. Shear/torsion fields use the Python `sectionproperties` service in Docker Compose; trusses and axial/bending-only states display without the service. **Member von Mises maxima · static** retains member-wide maxima; the [diagnostic page](reference/section-fem-lab.html) retains independent geometry and force input. [Viewing and limits](reference/EN/08-result-recovery.html#section-view-title).

## Getting started

Install Docker Desktop or Docker Engine with Compose. To use the current features, run from the cloned repository root:

```sh
docker compose up --build -d
```

Open [http://127.0.0.1:8080](http://127.0.0.1:8080). Stop with `docker compose down`.

In the preprocessor, **Load example…** beside **Examples** loads a complete model: two simple examples for each static, modal and transient analysis in 2D and 3D (12 total). Select one and click **Run analysis**. Examples also load when the editor is opened offline; edit the model and save it with **Download TXT**.

For the prebuilt release:

```sh
docker compose -f compose.release.yaml up -d
```

[compose.release.yaml](compose.release.yaml) pins all three services (solver, section solver and web interface) to **0.8.6**. Stop with `docker compose -f compose.release.yaml down`.

## Release verification

With PowerShell 7+ and Docker running, execute
`pwsh -NoProfile -File ./scripts/testing/test_release.ps1`. It builds the images,
runs the complete acceptance suite and writes
`output/release-runs/<run-id>/report.md`. [Verification guide](reference/EN/project-operations.html#checks).

[Script catalog (RU)](scripts/README.md) · [Test layout and registration (RU)](tests/README.md).

## Documentation and labs

- [Handbook](reference/EN/index.html) · [3D quick start](reference/EN/3d-quick-start.html) · [Example model](examples/cases/CaseSpatialQuickStart.txt).
- [Setup, configuration and checks](reference/EN/project-operations.html) · [Changelog (Russian)](CHANGELOG.md).
- [NAFEMS Challenge 5 benchmark](reference/EN/06a-nafems-challenge-5.html).

For the labs, install GNU Octave locally, open the repository folder, run `setup` and follow the handbook. Labs assume interactive work with code and matrices; Docker primarily serves the web interface.
