# Служебные команды

`scripts` содержит команды для сборки, расчётов, подготовки данных и запуска проверок. Сами тесты находятся в [tests](../tests/README.md).

| Каталог | Назначение |
|---|---|
| `build/` | Сборка ресурсов приложения и генерация каталога опор |
| `runtime/` | Неинтерактивный запуск Octave из backend |
| `examples/` | Экспорт примеров для интерфейса и справочника |
| `benchmarks/nafems/` | Генерация и проверка эталонной задачи NAFEMS Challenge 5 |
| `testing/` | Запуск наборов проверок, подготовка данных и релизный манифест |
| `release/` | Сборка и публикация версионированных Docker-образов |

Команды ниже выполняются из корня репозитория. PowerShell-сценарии запускаются через PowerShell 7 (`pwsh`); расширение `.ps1` обозначает среду исполнения, а не отдельную категорию Windows. Численные Python-проверки сечений выполняются в Linux-контейнерах.

| Задача | Команда | Среда | Docker | Результат |
|---|---|---|---|---|
| Сборка пространственного модуля | `npm run build:spatial` | Node.js, зависимости `npm ci` | Нет | `shared/vendor/` |
| Проверка актуальности сборки | `npm run check:spatial` | Node.js, зависимости `npm ci` | Нет | Вывод в терминал |
| Проверка каталога опор | `node scripts/build/generate-support-catalog.cjs --check` | Node.js | Нет | Вывод в терминал |
| Все численные проверки | `setup; run_octave_tests` | GNU Octave | Нет | `output/` |
| Экспорт примеров справочника | `setup; export_postprocessor_examples` | GNU Octave | Нет | `output/postprocessor-examples/` |
| Экспорт примеров MPC | `setup; export_mpc_examples` | GNU Octave | Нет | `output/` |
| Эталон NAFEMS | `setup; run_nafems_challenge5('quick')` (или `'full'`) | GNU Octave | Нет | `output/nafems-challenge-5/` |
| План релизных проверок | `pwsh -NoProfile -File scripts/testing/test_release.ps1 -PlanOnly` | PowerShell 7 | Нет | Проверка манифеста и план в терминале |
| Проверка механизма запуска | `pwsh -NoProfile -File scripts/testing/check_release_runner.ps1` | PowerShell 7 | Нет | `output/release-runner-checks/` |
| Полная приёмка | `pwsh -NoProfile -File scripts/testing/test_release.ps1` | PowerShell 7 | Да | `output/release-runs/<run-id>/` |
| Проверки восстановления результатов | `pwsh -NoProfile -File scripts/testing/test_recovery_pipeline.ps1` | PowerShell 7, готовые тестовые образы | Да | `output/recovery-runs/<run-id>/` |
| Проверки решателя сечений | `pwsh -NoProfile -File scripts/testing/run_section_solver_acceptance.ps1` | PowerShell 7 | Да | `output/section-solver/<run-id>/` |
| Проверка работающего стека | `pwsh -NoProfile -File scripts/testing/compose_smoke_test.ps1` | PowerShell 7, запущенный Compose-стек | Да | Вывод в терминал |
| Справка по сборке релиза | `python scripts/release/release_images.py --help` | Python 3 | Нет для справки; да для сборки | Вывод в терминал |

`runtime/run_solver_job.m` вызывается backend с путями входного и выходного файлов и названием результата. Путь переопределяется переменной `MKEF_RUNNER_SCRIPT`.

Полная приёмка собирает собственные образы, использует свежие данные и отдельный Compose-проект, сохраняет отчёты и удаляет свои временные контейнеры и образы. Публикация не выполняется. Подробности и выборочные проверки описаны в [руководстве разработчика](../DEVGUIDE.md).

Старые пути файлов в корне `scripts` удалены. Используйте новые пути; имена функций Octave и npm-команд сохранены.
