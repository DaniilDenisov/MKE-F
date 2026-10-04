'use strict';

    const stages = [
      {
        short: 'Модель',
        title: 'Модель и степени свободы',
        description: 'У каждого узла одна осевая степень свободы. Нумерация определяет адреса коэффициентов в глобальной матрице.',
        equation: 'K = 0'
      },
      {
        short: 'Элемент',
        title: 'Элементная матрица',
        description: 'Выберите любой стержень: его матрица зависит от k = EA/L, а локальные номера отображаются в две глобальные СС.',
        equation: 'k⁽ᵉ⁾ = (EA/L) · [1 −1; −1 1]'
      },
      {
        short: 'Вклад e₁',
        title: 'Сборка элемента e₁',
        description: 'Карта [u₁, u₂] переносит четыре коэффициента первой элементной матрицы в K.',
        equation: 'K([1,2],[1,2]) += k⁽¹⁾'
      },
      {
        short: 'Вклад e₂',
        title: 'Сборка элемента e₂',
        description: 'Вторая матрица попадает в строки и столбцы [u₂, u₃]. В K₂₂ два вклада складываются.',
        equation: 'K([2,3],[2,3]) += k⁽²⁾'
      },
      {
        short: 'Вклад e₃',
        title: 'Сборка элемента e₃',
        description: 'После третьего вклада глобальная матрица полностью собрана и имеет трёхдиагональный портрет.',
        equation: 'K([3,4],[3,4]) += k⁽³⁾'
      },
      {
        short: 'Заделка',
        title: 'Заделка u₁ = 0',
        description: 'Закреплённая степень свободы исключается. Решается свободный блок Kff для [u₂, u₃, u₄].',
        equation: 'Kff · uf = Ff'
      },
      {
        short: 'Решение',
        title: 'Перемещения и реакция',
        description: 'Свободные перемещения возвращаются в полный вектор, а реакция вычисляется по исходной системе: R = Ku − F.',
        equation: 'R = K · u − F'
      }
    ];

    const state = {
      stage: 0,
      selectedElement: 0,
      stiffness: [10, 20, 15],
      load: 30,
      view: 'dense',
      playing: false,
      timer: null
    };

    const dom = {
      stageStrip: document.getElementById('stageStrip'),
      stageTitle: document.getElementById('stageTitle'),
      stageDescription: document.getElementById('stageDescription'),
      stageBadge: document.getElementById('stageBadge'),
      prevButton: document.getElementById('prevButton'),
      nextButton: document.getElementById('nextButton'),
      playButton: document.getElementById('playButton'),
      resetButton: document.getElementById('resetButton'),
      speedRange: document.getElementById('speedRange'),
      elementPicker: document.getElementById('elementPicker'),
      localTitle: document.getElementById('localTitle'),
      localStiffness: document.getElementById('localStiffness'),
      dofMapChip: document.getElementById('dofMapChip'),
      localMatrixBody: document.getElementById('localMatrixBody'),
      globalMatrix: document.getElementById('globalMatrix'),
      matrixTitle: document.getElementById('matrixTitle'),
      equationChip: document.getElementById('equationChip'),
      denseView: document.getElementById('denseView'),
      sparseView: document.getElementById('sparseView'),
      tripletList: document.getElementById('tripletList'),
      csrTitle: document.getElementById('csrTitle'),
      csrValues: document.getElementById('csrValues'),
      csrColumns: document.getElementById('csrColumns'),
      csrRowPtr: document.getElementById('csrRowPtr'),
      storageComparison: document.getElementById('storageComparison'),
      operationLine: document.getElementById('operationLine'),
      resultGrid: document.getElementById('resultGrid'),
      footnote: document.getElementById('footnote'),
      deformedShape: document.getElementById('deformedShape'),
      reactionArrow: document.getElementById('reactionArrow'),
      reactionLabel: document.getElementById('reactionLabel'),
      loadLabel: document.getElementById('loadLabel')
    };

    function formatNumber(value, digits = 3) {
      if (Math.abs(value) < 1e-10) return '0';
      const rounded = Number(value.toFixed(digits));
      return String(rounded).replace('-', '−');
    }

    function zeroMatrix(size) {
      return Array.from({ length: size }, () => Array(size).fill(0));
    }

    function assembledElementCount() {
      return Math.max(0, Math.min(3, state.stage - 1));
    }

    function fullMatrix(count = assembledElementCount()) {
      const matrix = zeroMatrix(4);
      for (let element = 0; element < count; element += 1) {
        const k = state.stiffness[element];
        matrix[element][element] += k;
        matrix[element][element + 1] -= k;
        matrix[element + 1][element] -= k;
        matrix[element + 1][element + 1] += k;
      }
      return matrix;
    }

    function solveLinearSystem(matrix, vector) {
      const n = vector.length;
      const augmented = matrix.map((row, index) => [...row, vector[index]]);
      for (let pivot = 0; pivot < n; pivot += 1) {
        let maxRow = pivot;
        for (let row = pivot + 1; row < n; row += 1) {
          if (Math.abs(augmented[row][pivot]) > Math.abs(augmented[maxRow][pivot])) maxRow = row;
        }
        [augmented[pivot], augmented[maxRow]] = [augmented[maxRow], augmented[pivot]];
        const divisor = augmented[pivot][pivot];
        if (Math.abs(divisor) < 1e-12) throw new Error('Матрица вырождена');
        for (let column = pivot; column <= n; column += 1) augmented[pivot][column] /= divisor;
        for (let row = 0; row < n; row += 1) {
          if (row === pivot) continue;
          const factor = augmented[row][pivot];
          for (let column = pivot; column <= n; column += 1) {
            augmented[row][column] -= factor * augmented[pivot][column];
          }
        }
      }
      return augmented.map(row => row[n]);
    }

    function solution() {
      const K = fullMatrix(3);
      const Kff = K.slice(1).map(row => row.slice(1));
      const uf = solveLinearSystem(Kff, [0, 0, state.load]);
      const u = [0, ...uf];
      const force = [0, 0, 0, state.load];
      const reaction = K.map((row, i) => row.reduce((sum, value, j) => sum + value * u[j], 0) - force[i]);
      return { K, Kff, uf, u, force, reaction };
    }

    function displayedMatrix() {
      if (state.stage >= 5) return solution().Kff;
      return fullMatrix();
    }

    function displayedLabels() {
      return state.stage >= 5 ? ['u₂', 'u₃', 'u₄'] : ['u₁', 'u₂', 'u₃', 'u₄'];
    }

    function currentContributionCells() {
      if (state.stage < 2 || state.stage > 4) return [];
      const element = state.stage - 2;
      return [
        `${element}-${element}`,
        `${element}-${element + 1}`,
        `${element + 1}-${element}`,
        `${element + 1}-${element + 1}`
      ];
    }

    function cellExpression(row, column, count) {
      const terms = [];
      for (let element = 0; element < count; element += 1) {
        if ((row === element || row === element + 1) && (column === element || column === element + 1)) {
          const positive = row === column;
          terms.push(`${positive ? '' : '−'}k${element + 1}`);
        }
      }
      return terms.length ? terms.join(' + ').replace('+ −', '− ') : '0';
    }

    function renderStageStrip() {
      dom.stageStrip.innerHTML = stages.map((stage, index) => `
        <button class="stage-button${index < state.stage ? ' is-complete' : ''}"
                type="button"
                data-stage="${index}"
                ${index === state.stage ? 'aria-current="step"' : ''}>
          <span>${index < state.stage ? '✓ завершено' : `шаг ${index + 1}`}</span>
          ${stage.short}
        </button>
      `).join('');
    }

    function renderElementPicker() {
      dom.elementPicker.innerHTML = state.stiffness.map((value, index) => `
        <button class="element-button" type="button" data-element-button="${index}"
                aria-pressed="${index === state.selectedElement}">
          e${index + 1} · k${index + 1} = ${formatNumber(value)}
        </button>
      `).join('');
    }

    function renderLocalMatrix() {
      const element = state.selectedElement;
      const k = state.stiffness[element];
      dom.localTitle.textContent = `Элемент e${element + 1}`;
      dom.localStiffness.textContent = `k${element + 1} = E${element + 1}A${element + 1}/L${element + 1} = ${formatNumber(k)} кН/мм`;
      dom.dofMapChip.textContent = `глобальные [u${element + 1}, u${element + 2}]`;
      dom.localMatrixBody.innerHTML = `
        <tr><td>${formatNumber(k)}</td><td>${formatNumber(-k)}</td></tr>
        <tr><td>${formatNumber(-k)}</td><td>${formatNumber(k)}</td></tr>
      `;
    }

    function renderModel() {
      const count = assembledElementCount();
      for (let index = 0; index < 3; index += 1) {
        const group = document.getElementById(`memberGroup${index + 1}`);
        const line = group.querySelector('.member-base');
        line.className.baseVal = 'member-base';
        if (index < count) line.classList.add('member-assembled');
        if (index === state.selectedElement) line.classList.add('member-active');
        if (state.stage >= 2 && state.stage <= 4 && index === state.stage - 2) line.classList.add('member-current');
      }

      dom.loadLabel.textContent = `P = ${formatNumber(state.load)} кН`;
      const solved = state.stage === 6;
      dom.deformedShape.toggleAttribute('hidden', !solved);
      dom.reactionArrow.toggleAttribute('hidden', !solved);
      dom.reactionLabel.textContent = `R₁ = −${formatNumber(state.load)} кН`;

      if (solved) {
        const values = solution().u;
        const maximum = Math.max(...values.map(Math.abs), 1);
        const points = [100, 300, 500, 700].map((x, index) => `${x + (values[index] / maximum) * 70},105`).join(' ');
        dom.deformedShape.querySelector('polyline').setAttribute('points', points);
        dom.deformedShape.querySelectorAll('circle').forEach((circle, index) => {
          circle.setAttribute('cx', String([100, 300, 500, 700][index] + (values[index] / maximum) * 70));
        });
      }
    }

    function renderDenseMatrix() {
      const matrix = displayedMatrix();
      const labels = displayedLabels();
      const contributionCells = new Set(currentContributionCells());
      const count = assembledElementCount();
      let html = '<thead><tr><th></th>' + labels.map(label => `<th scope="col">${label}</th>`).join('') + '</tr></thead><tbody>';

      matrix.forEach((row, rowIndex) => {
        html += `<tr><th scope="row">${labels[rowIndex]}</th>`;
        row.forEach((value, columnIndex) => {
          const isReduced = state.stage >= 5;
          const sourceRow = isReduced ? rowIndex + 1 : rowIndex;
          const sourceColumn = isReduced ? columnIndex + 1 : columnIndex;
          const key = `${sourceRow}-${sourceColumn}`;
          const classes = [Math.abs(value) > 1e-10 ? 'nonzero' : 'zero'];
          if (contributionCells.has(key)) classes.push('current-cell');
          const expression = state.stage >= 5 ? cellExpression(sourceRow, sourceColumn, 3) : cellExpression(rowIndex, columnIndex, count);
          html += `<td class="${classes.join(' ')}" title="${labels[rowIndex]}, ${labels[columnIndex]}: ${expression}">
            ${formatNumber(value)}
            ${Math.abs(value) > 1e-10 ? `<span class="cell-note">${expression}</span>` : ''}
          </td>`;
        });
        html += '</tr>';
      });
      html += '</tbody>';
      dom.globalMatrix.innerHTML = html;
      dom.globalMatrix.setAttribute('aria-label', state.stage >= 5 ? 'Редуцированная матрица жёсткости' : 'Глобальная матрица жёсткости');
    }

    function rawTriplets() {
      const count = state.stage >= 5 ? 3 : assembledElementCount();
      const triplets = [];
      for (let element = 0; element < count; element += 1) {
        const k = state.stiffness[element];
        triplets.push(
          { row: element + 1, column: element + 1, value: k, element },
          { row: element + 1, column: element + 2, value: -k, element },
          { row: element + 2, column: element + 1, value: -k, element },
          { row: element + 2, column: element + 2, value: k, element }
        );
      }
      if (state.stage >= 5) {
        return triplets
          .filter(item => item.row !== 1 && item.column !== 1)
          .map(item => ({ ...item, row: item.row - 1, column: item.column - 1 }));
      }
      return triplets;
    }

    function matrixToCsr(matrix) {
      const values = [];
      const columns = [];
      const rowPtr = [1];
      matrix.forEach(row => {
        row.forEach((value, column) => {
          if (Math.abs(value) > 1e-10) {
            values.push(value);
            columns.push(column + 1);
          }
        });
        rowPtr.push(values.length + 1);
      });
      return { values, columns, rowPtr };
    }

    function renderSparseMatrix() {
      const triplets = rawTriplets();
      const currentElement = state.stage >= 2 && state.stage <= 4 ? state.stage - 2 : -1;
      dom.tripletList.innerHTML = triplets.length
        ? triplets.map(item => `<span class="triplet${item.element === currentElement ? ' is-new' : ''}">(${item.row}, ${item.column}, ${formatNumber(item.value)})</span>`).join('')
        : '<span class="triplet">пока пусто</span>';

      const matrix = displayedMatrix();
      const csr = matrixToCsr(matrix);
      dom.csrTitle.textContent = `После sparse() · ${matrix.length} строк · учебное чтение по строкам`;
      dom.csrValues.textContent = `[${csr.values.map(value => formatNumber(value)).join(', ')}]`;
      dom.csrColumns.textContent = `[${csr.columns.join(', ')}]`;
      dom.csrRowPtr.textContent = `[${csr.rowPtr.join(', ')}]`;
      const denseNumbers = matrix.length * matrix.length;
      dom.storageComparison.textContent = `Плотно: ${denseNumbers} чисел · разреженно: ${csr.values.length} ненулевых значений; columns и rowPtr ниже показаны только как учебное построчное представление`;
    }

    function operationText() {
      const k = state.stiffness;
      const P = state.load;
      if (state.stage === 0) return 'Пока K заполнена нулями. Карта элемента определит, в какие четыре ячейки попадёт его вклад.';
      if (state.stage === 1) return `Для выбранного e${state.selectedElement + 1}: k = EA/L = ${formatNumber(k[state.selectedElement])} кН/мм; локальные коэффициенты ещё не добавлены в K.`;
      if (state.stage === 2) return `K₁₁ += k₁; K₁₂ += −k₁; K₂₁ += −k₁; K₂₂ += k₁. Получено K₂₂ = ${formatNumber(k[0])}.`;
      if (state.stage === 3) return `Повторяющийся адрес суммируется: K₂₂ = k₁ + k₂ = ${formatNumber(k[0])} + ${formatNumber(k[1])} = ${formatNumber(k[0] + k[1])}.`;
      if (state.stage === 4) return `K₃₃ = k₂ + k₃ = ${formatNumber(k[1])} + ${formatNumber(k[2])} = ${formatNumber(k[1] + k[2])}. Глобальная K собрана.`;
      if (state.stage === 5) return `u₁ = 0, поэтому строка и столбец 1 не входят в Kff. Ff = [0, 0, ${formatNumber(P)}]ᵀ кН; исходная K сохранена для реакций.`;
      const u = solution().u;
      return `u = [${u.map(value => formatNumber(value)).join(', ')}]ᵀ мм; R₁ = −${formatNumber(P)} кН; каждый элемент передаёт N = ${formatNumber(P)} кН.`;
    }

    function renderResults() {
      const count = assembledElementCount();
      if (state.stage === 6) {
        const result = solution();
        dom.resultGrid.innerHTML = `
          <div class="result-card">
            <span class="result-card__label">Перемещения u, мм</span>
            <span class="result-card__value">[${result.u.map(value => formatNumber(value)).join(', ')}]ᵀ</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">Реакция R₁</span>
            <span class="result-card__value">−${formatNumber(state.load)} кН</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">Продольные силы</span>
            <span class="result-card__value">N₁ = N₂ = N₃ = ${formatNumber(state.load)} кН</span>
          </div>
        `;
      } else {
        dom.resultGrid.innerHTML = `
          <div class="result-card">
            <span class="result-card__label">Собрано элементов</span>
            <span class="result-card__value">${count} / 3</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">Размер решаемой системы</span>
            <span class="result-card__value">${state.stage >= 5 ? '3 × 3 (из 4 × 4)' : '4 × 4'}</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">Состояние</span>
            <span class="result-card__value">${state.stage === 0 ? 'K не собрана' : state.stage === 1 ? 'Матрица элемента' : state.stage < 4 ? 'Сборка идёт' : state.stage === 4 ? 'K собрана' : 'Kff готова'}</span>
          </div>
        `;
      }

      dom.footnote.textContent = state.stage >= 5
        ? 'Редукция показана как исключение закреплённой СС. Полная K нужна для восстановления реакции R = Ku − F.'
        : 'В плотном виде видны все ячейки, включая нули. В программе триплеты передаются библиотечной функции GNU Octave sparse(), которая суммирует совпадающие адреса и создаёт разреженную матрицу.';
    }

    function renderViewMode() {
      const dense = state.view === 'dense';
      dom.denseView.hidden = !dense;
      dom.sparseView.hidden = dense;
      document.querySelectorAll('[data-view]').forEach(button => {
        button.setAttribute('aria-pressed', String(button.dataset.view === state.view));
      });
    }

    function render() {
      const stage = stages[state.stage];
      dom.stageTitle.textContent = stage.title;
      dom.stageDescription.textContent = stage.description;
      dom.stageBadge.textContent = `Шаг ${state.stage + 1} из ${stages.length}`;
      dom.equationChip.textContent = stage.equation;
      dom.matrixTitle.textContent = state.stage >= 5 ? 'Редуцированная матрица Kff' : 'Глобальная матрица жёсткости K';
      dom.operationLine.textContent = operationText();
      dom.prevButton.disabled = state.stage === 0;
      dom.nextButton.disabled = state.stage === stages.length - 1;
      dom.playButton.textContent = state.playing ? '❚❚ Пауза' : '▶ Запустить';

      renderStageStrip();
      renderElementPicker();
      renderLocalMatrix();
      renderModel();
      renderDenseMatrix();
      renderSparseMatrix();
      renderResults();
      renderViewMode();
    }

    function selectStage(index) {
      state.stage = Math.max(0, Math.min(stages.length - 1, index));
      if (state.stage >= 2 && state.stage <= 4) state.selectedElement = state.stage - 2;
      if (state.stage === stages.length - 1 && state.playing) stopPlayback();
      render();
    }

    function stopPlayback() {
      state.playing = false;
      window.clearInterval(state.timer);
      state.timer = null;
    }

    function restartPlaybackTimer() {
      window.clearInterval(state.timer);
      const delays = { 1: 2600, 2: 1700, 3: 1050 };
      state.timer = window.setInterval(() => {
        if (state.stage >= stages.length - 1) {
          stopPlayback();
          render();
          return;
        }
        selectStage(state.stage + 1);
      }, delays[dom.speedRange.value]);
    }

    function togglePlayback() {
      if (state.playing) {
        stopPlayback();
        render();
        return;
      }
      if (state.stage === stages.length - 1) state.stage = 0;
      state.playing = true;
      render();
      restartPlaybackTimer();
    }

    function selectElement(index) {
      state.selectedElement = index;
      if (state.stage === 0) state.stage = 1;
      render();
    }

    dom.prevButton.addEventListener('click', () => selectStage(state.stage - 1));
    dom.nextButton.addEventListener('click', () => selectStage(state.stage + 1));
    dom.playButton.addEventListener('click', togglePlayback);
    dom.resetButton.addEventListener('click', () => {
      stopPlayback();
      state.stage = 0;
      state.selectedElement = 0;
      render();
    });

    dom.speedRange.addEventListener('input', () => {
      if (state.playing) restartPlaybackTimer();
    });

    document.addEventListener('click', event => {
      const stageButton = event.target.closest('[data-stage]');
      if (stageButton) {
        stopPlayback();
        selectStage(Number(stageButton.dataset.stage));
        return;
      }

      const viewButton = event.target.closest('[data-view]');
      if (viewButton) {
        state.view = viewButton.dataset.view;
        renderViewMode();
        return;
      }

      const elementButton = event.target.closest('[data-element-button]');
      if (elementButton) {
        selectElement(Number(elementButton.dataset.elementButton));
        return;
      }

      const elementGroup = event.target.closest('[data-element]');
      if (elementGroup) selectElement(Number(elementGroup.dataset.element));
    });

    document.querySelectorAll('[data-element]').forEach(group => {
      group.addEventListener('keydown', event => {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          selectElement(Number(group.dataset.element));
        }
      });
    });

    document.querySelectorAll('[data-stiffness]').forEach(input => {
      input.addEventListener('input', () => {
        const index = Number(input.dataset.stiffness);
        state.stiffness[index] = Number(input.value);
        document.getElementById(`stiffnessValue${index + 1}`).textContent = `${formatNumber(state.stiffness[index])} кН/мм`;
        render();
      });
    });

    document.getElementById('loadRange').addEventListener('input', event => {
      state.load = Number(event.target.value);
      document.getElementById('loadValue').textContent = `${formatNumber(state.load)} кН`;
      render();
    });

    render();
