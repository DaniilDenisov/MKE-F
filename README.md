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
staticResult = problem.RunStatic();
```

Методы `RunStatic`, `RunModal` и `RunTransient` возвращают структуры с явно
названными результатами и не изменяют собранные матрицы `K` и `M` или старое
совместимое поле `F`. Численное ядро можно вызывать без печати и графиков:

```octave
model = problem.GetAnalysisModel();
staticResult = solveStatic(model);
modalResult = solveModal(model);
dynamicResult = solveTransient(model, ...
    struct('timeStep', 1e-4, 'duration', 1e-2));
```

Опоры преобразуются в списки `fixedDOFs` и `freeDOFs`. Решатели работают с
редуцированными матрицами свободных степеней свободы, а возвращаемые полные
векторы и формы колебаний содержат точные нули на закреплённых степенях свободы.
Статический результат также содержит реакции и невязку общего равновесия
`equilibriumResidual = [Fx; Fy; Mz]`.

## Узловые нагрузки

Для плоской фермы узловая нагрузка имеет компоненты `[Fx, Fy]`. Для плоской
рамы компоненты `[Fx, Fy, Mz]` соответствуют степеням свободы `[ux, uy, θz]`;
силы `Fz` в двумерной модели нет.

- `bcforce_stat`, тип `10`: статическая нагрузка в `RunStatic`; в
  `RunTransient` сохраняет старое значение однократного прямоугольного импульса
  на первом шаге, с дискретным импульсом `F0*dt`.
- `bcforce_harm`, тип `11`: гармоническая нагрузка `F0*sin(2*pi*f*t)`.
- `bcforce_pulse`, тип `12`: явно заданный однократный импульс на первом шаге.
- `bcforce_step`, тип `13`: постоянная ступенчатая нагрузка на всех шагах.

Несколько строк нагрузки в одном блоке суммируются, в том числе если они
относятся к одному узлу или одной степени свободы.

## Переходный расчёт

Все возвращаемые истории используют одну временную сетку: столбец 1
соответствует `t=0`, столбец 2 — `t=dt`, последний столбец — фактической
длительности расчёта. Это относится к `time`, `loadHistory`, `displacements`,
`velocities` и `accelerations`.

Начальные перемещения и скорости можно передать численному ядру:

```octave
transientOptions = struct(...
    'timeStep', 1e-4, ...
    'duration', 1e-2, ...
    'initialDisplacement', u0, ...
    'initialVelocity', v0);
result = solveTransient(model, transientOptions);
```

Или через совместимый фасад:

```octave
initial = struct('initialDisplacement', u0, 'initialVelocity', v0);
result = problem.RunTransient(1e-4, 1e-2, 2, 1, initial);
```

Начальное ускорение вычисляется из равновесия. Результат также содержит
`reactions`, `equilibriumResidual`, `spectrumFrequencyHz` и
`displacementAmplitudeSpectrum`. Используется метод Ньюмарка со средней
акселерацией (`beta=1/4`, `gamma=1/2`) по формулам (22)-(23) из конспекта
[H. P. Gavin, Numerical Integration in Structural Dynamics](https://people.duke.edu/~hpgavin/StructuralDynamics/NumericalIntegration.pdf).
