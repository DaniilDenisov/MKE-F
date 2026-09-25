# МКЭ-Ф

Учебная программа на MATLAB/GNU Octave для расчета плоских ферм и рам
методом конечных элементов.

https://конструкторский.рф/

## GNU Octave

Минимальная целевая версия — GNU Octave 8.4. Локально программа проверяется
на GNU Octave 11.3.0, а более старая версия проверяется в CI на Ubuntu 24.04.

Запуск переносимого smoke-теста без графического интерфейса:

```sh
octave --no-gui --quiet --eval "addpath(pwd); run_octave_smoke_tests;"
```

Для запуска задачи без построения графиков и диагностической печати:

```matlab
options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem('Case1ElementBeam.txt', options);
problem.RunStatic();
```
