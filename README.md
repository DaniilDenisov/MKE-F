# МКЭ-Ф

Учебная программа на GNU Octave для расчета плоских ферм и рам
методом конечных элементов.

https://конструкторский.рф/

## GNU Octave

Минимальная целевая версия — GNU Octave 8.4. Локально программа проверяется
на GNU Octave 11.3.0, а более старая версия проверяется в CI на Ubuntu 24.04.

GNU Octave является единственной поддерживаемой средой выполнения. Совместимость
с MATLAB не тестируется и не гарантируется.

Запуск полного набора тестов без графического интерфейса:

```sh
octave --no-gui --quiet --eval "addpath(pwd); run_octave_tests;"
```

Запуск только быстрого smoke-теста:

```sh
octave --no-gui --quiet --eval "addpath(pwd); run_octave_smoke_tests;"
```

Для запуска задачи без построения графиков и диагностической печати:

```octave
options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem('Case1ElementBeam.txt', options);
problem.RunStatic();
```
