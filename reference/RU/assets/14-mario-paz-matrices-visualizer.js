'use strict';

    const stages = [
      { short: 'Схема', title: 'Расчётная схема и степени свободы', description: 'Три рамных узла дают девять глобальных степеней свободы. Наружные узлы закреплены, средний узел свободен.', equation: '9 СС = 3 узла × 3' },
      { short: 'Локальные', title: 'Локальные матрицы элемента', description: 'Оба элемента имеют одинаковые A, E, ρ, I и L. Поэтому их локальные матрицы жёсткости одинаковы между собой; то же справедливо для локальных матриц масс.', equation: 'kˡ, mˡ · 6 × 6' },
      { short: 'Поворот e₁', title: 'Поворот наклонного элемента e₁', description: 'Для c = s = 1/√2 локальные матрицы переводятся в глобальные оси одинаковым преобразованием T₁.', equation: 'k₁ = T₁ᵀkˡT₁' },
      { short: 'Элемент e₂', title: 'Горизонтальный элемент e₂', description: 'Для второго элемента c = 1 и s = 0, поэтому T₂ = I: локальные и глобальные компоненты совпадают.', equation: 'k₂ = kˡ, m₂ = mˡ' },
      { short: 'Вклад e₁', title: 'Первый вклад в систему 9 × 9', description: 'Матрица e₁ добавляется по карте [1, 2, 3, 4, 5, 6]. Выделенные ячейки принадлежат текущему элементу.', equation: 'K(d₁,d₁) += k₁' },
      { short: 'Вклад e₂', title: 'Сборка двух элементов', description: 'Матрица e₂ добавляется по карте [4, 5, 6, 7, 8, 9]. В общем узле 2 вклады двух элементов суммируются.', equation: 'K(d₂,d₂) += k₂' },
      { short: 'Заделки', title: 'Редукция до свободного узла', description: 'Строки и столбцы закреплённых СС [1, 2, 3, 7, 8, 9] исключаются. Остаются Kff и Mff размера 3 × 3.', equation: 'free = [4, 5, 6]' },
      { short: 'Динамика', title: 'Система переходного анализа', description: 'Редуцированные Kff и Mff образуют уравнение движения. Метод Ньюмарка решает систему с постоянной эффективной матрицей.', equation: 'Mff·ü + Kff·u = F(t)' }
    ];

    const state = {
      stage: 0,
      matrix: 'K',
      selectedElement: 0,
      selectedCell: null,
      timeStep: 0.01,
      playing: false,
      timer: null
    };

    const properties = { A: 6, E: 1e7, rho: 0.7, I: 100, L: 100 };
    const free = [3, 4, 5];
    const fixed = [0, 1, 2, 6, 7, 8];
    const elements = [
      { name: 'e₁', c: Math.SQRT1_2, s: Math.SQRT1_2, map: [0, 1, 2, 3, 4, 5], nodes: [1, 2] },
      { name: 'e₂', c: 1, s: 0, map: [3, 4, 5, 6, 7, 8], nodes: [2, 3] }
    ];

    const fullLabels = ['uₓ₁', 'uᵧ₁', 'θ₁', 'uₓ₂', 'uᵧ₂', 'θ₂', 'uₓ₃', 'uᵧ₃', 'θ₃'];
    const localLabels = ['u₁', 'v₁', 'θ₁', 'u₂', 'v₂', 'θ₂'];
    const reducedLabels = ['uₓ₂', 'uᵧ₂', 'θ₂'];

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
      dofGrid: document.getElementById('dofGrid'),
      formulaTitle: document.getElementById('formulaTitle'),
      formulaLines: document.getElementById('formulaLines'),
      matrixTitle: document.getElementById('matrixTitle'),
      matrixSubtitle: document.getElementById('matrixSubtitle'),
      equationChip: document.getElementById('equationChip'),
      matrixTable: document.getElementById('matrixTable'),
      cellDetail: document.getElementById('cellDetail'),
      summaryGrid: document.getElementById('summaryGrid'),
      newmarkPanel: document.getElementById('newmarkPanel'),
      timeStepRange: document.getElementById('timeStepRange'),
      timeStepValue: document.getElementById('timeStepValue'),
      effectiveFormula: document.getElementById('effectiveFormula'),
      dynamicDofs: document.getElementById('dynamicDofs')
    };

    function zeros(rows, columns = rows) {
      return Array.from({ length: rows }, () => Array(columns).fill(0));
    }

    function transpose(matrix) {
      return matrix[0].map((_, column) => matrix.map(row => row[column]));
    }

    function multiply(left, right) {
      const result = zeros(left.length, right[0].length);
      for (let row = 0; row < left.length; row += 1) {
        for (let column = 0; column < right[0].length; column += 1) {
          for (let inner = 0; inner < right.length; inner += 1) result[row][column] += left[row][inner] * right[inner][column];
        }
      }
      return result;
    }

    function add(left, right) {
      return left.map((row, i) => row.map((value, j) => value + right[i][j]));
    }

    function scaleMatrix(matrix, factor) {
      return matrix.map(row => row.map(value => value * factor));
    }

    function submatrix(matrix, indices) {
      return indices.map(row => indices.map(column => matrix[row][column]));
    }

    function localMatrices() {
      const { A, E, rho, I, L } = properties;
      const axial = E * A / L;
      const b12 = 12 * E * I / L ** 3;
      const b6 = 6 * E * I / L ** 2;
      const b4 = 4 * E * I / L;
      const b2 = 2 * E * I / L;
      const K = [
        [axial, 0, 0, -axial, 0, 0],
        [0, b12, b6, 0, -b12, b6],
        [0, b6, b4, 0, -b6, b2],
        [-axial, 0, 0, axial, 0, 0],
        [0, -b12, -b6, 0, b12, -b6],
        [0, b6, b2, 0, -b6, b4]
      ];
      const factor = rho * A * L / 420;
      const M = scaleMatrix([
        [140, 0, 0, 70, 0, 0],
        [0, 156, 22 * L, 0, 54, -13 * L],
        [0, 22 * L, 4 * L ** 2, 0, 13 * L, -3 * L ** 2],
        [70, 0, 0, 140, 0, 0],
        [0, 54, 13 * L, 0, 156, -22 * L],
        [0, -13 * L, -3 * L ** 2, 0, -22 * L, 4 * L ** 2]
      ], factor);
      return { K, M };
    }

    function transformation(c, s) {
      return [
        [c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]
      ];
    }

    function elementMatrix(index, kind) {
      const local = localMatrices()[kind];
      const T = transformation(elements[index].c, elements[index].s);
      return multiply(multiply(transpose(T), local), T);
    }

    function assembledMatrix(kind, count = 2) {
      const global = zeros(9);
      for (let element = 0; element < count; element += 1) {
        const matrix = elementMatrix(element, kind);
        const map = elements[element].map;
        for (let i = 0; i < 6; i += 1) {
          for (let j = 0; j < 6; j += 1) global[map[i]][map[j]] += matrix[i][j];
        }
      }
      return global;
    }

    function reducedMatrices() {
      return {
        K: submatrix(assembledMatrix('K'), free),
        M: submatrix(assembledMatrix('M'), free)
      };
    }

    function effectiveMatrix() {
      const reduced = reducedMatrices();
      const a0 = 4 / state.timeStep ** 2;
      return add(reduced.K, scaleMatrix(reduced.M, a0));
    }

    function displayContext() {
      let kind = state.matrix;
      if (kind === 'Keff' && state.stage < 7) kind = 'K';
      const baseKind = kind === 'Keff' ? 'Keff' : kind;

      if (state.stage === 0) return { matrix: zeros(9), labels: fullLabels, kind: baseKind, title: kindTitle(kind, 'Глобальная'), subtitle: 'нулевая заготовка · 9 × 9', scale: matrixScale(kind) };
      if (state.stage === 1) return { matrix: localMatrices()[kind], labels: localLabels, kind, title: kindTitle(kind, 'Локальная'), subtitle: 'порядок [u₁, v₁, θ₁, u₂, v₂, θ₂] · 6 × 6', scale: matrixScale(kind) };
      if (state.stage === 2 || state.stage === 3) {
        const element = state.stage - 2;
        const labels = elements[element].map.map(index => fullLabels[index]);
        return { matrix: elementMatrix(element, kind), labels, kind, title: `${kindTitle(kind, 'Матрица')} элемента e${element + 1}`, subtitle: `${element === 0 ? 'T₁ᵀ · local · T₁' : 'T₂ = I'} · 6 × 6`, scale: matrixScale(kind) };
      }
      if (state.stage === 4 || state.stage === 5) {
        const count = state.stage - 3;
        return { matrix: assembledMatrix(kind, count), labels: fullLabels, kind, title: kindTitle(kind, 'Глобальная'), subtitle: `${count} из 2 элементов собрано · 9 × 9`, scale: matrixScale(kind) };
      }
      if (state.stage === 6) {
        return {
          matrix: assembledMatrix(kind),
          labels: fullLabels,
          kind,
          title: kindTitle(kind, 'Редуцированная'),
          subtitle: 'исходная 9 × 9 · серым погашены закреплённые строки и столбцы · активный блок 3 × 3',
          scale: matrixScale(kind),
          reductionPreview: true
        };
      }
      if (kind === 'Keff') return { matrix: effectiveMatrix(), labels: reducedLabels, kind, title: 'Эффективная матрица Kₑff', subtitle: `β = 1/4 · Δt = ${state.timeStep.toFixed(3)} s · 3 × 3`, scale: 1e6 };
      const reduced = reducedMatrices()[kind];
      return { matrix: reduced, labels: reducedLabels, kind, title: kindTitle(kind, 'Редуцированная'), subtitle: 'свободные СС [4, 5, 6] · 3 × 3', scale: matrixScale(kind) };
    }

    function kindTitle(kind, prefix) {
      if (kind === 'M') return `${prefix} матрица масс M`;
      if (kind === 'Keff') return 'Эффективная матрица Kₑff';
      return `${prefix} матрица жёсткости K`;
    }

    function matrixScale(kind) {
      return kind === 'M' ? 1 : 1e6;
    }

    function formatNumber(value, digits = 4) {
      if (Math.abs(value) < 5e-10) return '0';
      const rounded = Number(value.toFixed(digits));
      return String(rounded).replace('-', '−');
    }

    function formatPhysical(value, kind) {
      if (kind === 'M') return `${formatNumber(value, 4)}`;
      return `${formatNumber(value / 1e6, 6)} × 10⁶`;
    }

    function renderStageStrip() {
      dom.stageStrip.innerHTML = stages.map((stage, index) => `
        <button class="stage-button${index < state.stage ? ' is-complete' : ''}" type="button" data-stage="${index}" ${index === state.stage ? 'aria-current="step"' : ''}>
          <span>${index < state.stage ? '✓ завершено' : `шаг ${index + 1}`}</span>${stage.short}
        </button>
      `).join('');
    }

    function renderElementPicker() {
      dom.elementPicker.innerHTML = elements.map((element, index) => `
        <button class="element-button" type="button" data-element-button="${index}" aria-pressed="${state.selectedElement === index}">
          ${element.name} · карта [${element.map.map(value => value + 1).join(', ')}]
        </button>
      `).join('');
    }

    function renderModel() {
      for (let index = 0; index < 2; index += 1) {
        const line = document.querySelector(`#elementGroup${index + 1} .member`);
        line.className.baseVal = 'member';
        if ((state.stage === 4 && index === 0) || (state.stage === 5 && index < 2) || state.stage >= 6) line.classList.add('member-assembled');
        if (index === state.selectedElement) line.classList.add('member-selected');
        if ((state.stage === 2 && index === 0) || (state.stage === 3 && index === 1) || (state.stage === 4 && index === 0) || (state.stage === 5 && index === 1)) line.classList.add('member-current');
      }
      dom.dynamicDofs.toggleAttribute('hidden', state.stage < 6);
    }

    function renderDofs() {
      dom.dofGrid.innerHTML = fullLabels.map((label, index) => {
        const constrained = state.stage >= 6 && fixed.includes(index);
        const freeDof = state.stage >= 6 && free.includes(index);
        return `<span class="dof-chip${constrained ? ' fixed' : ''}${freeDof ? ' free' : ''}">${index + 1}: ${label}</span>`;
      }).join('');
    }

    function renderFormulaPanel() {
      if (state.matrix === 'M') {
        dom.formulaTitle.textContent = 'Согласованная матрица масс';
        dom.formulaLines.innerHTML = `
          <div class="formula-line">ρAL / 420 = 0.7 · 6 · 100 / 420 = 1</div>
          <div class="formula-line">mᵉ = (ρAL/420) · m̄(L)</div>
          <div class="formula-line">mᵍ = Tᵀ · mᵉ · T</div>`;
      } else if (state.matrix === 'Keff') {
        const a0 = 4 / state.timeStep ** 2;
        dom.formulaTitle.textContent = 'Эффективная жёсткость Ньюмарка';
        dom.formulaLines.innerHTML = `
          <div class="formula-line">β = 1/4, γ = 1/2</div>
          <div class="formula-line">a₀ = 1/(βΔt²) = ${formatNumber(a0, 3)}</div>
          <div class="formula-line">Keff = Kff + a₀ · Mff</div>`;
      } else {
        dom.formulaTitle.textContent = 'Коэффициенты локальной матрицы жёсткости';
        dom.formulaLines.innerHTML = `
          <div class="formula-line">EA/L = 600 000</div>
          <div class="formula-line">12EI/L³ = 12 000 · 6EI/L² = 600 000</div>
          <div class="formula-line">4EI/L = 40 000 000 · 2EI/L = 20 000 000</div>`;
      }
    }

    function currentGlobalCells() {
      if (state.stage !== 4 && state.stage !== 5) return new Set();
      const element = state.stage - 4;
      const matrix = elementMatrix(element, state.matrix === 'M' ? 'M' : 'K');
      const map = elements[element].map;
      const result = new Set();
      for (let i = 0; i < 6; i += 1) {
        for (let j = 0; j < 6; j += 1) {
          if (Math.abs(matrix[i][j]) > 1e-9) result.add(`${map[i]}-${map[j]}`);
        }
      }
      return result;
    }

    function renderMatrix() {
      const context = displayContext();
      const current = currentGlobalCells();
      const compact = context.matrix.length >= 9;
      dom.matrixTable.classList.toggle('compact', compact);
      dom.matrixTitle.textContent = context.title;
      dom.matrixSubtitle.textContent = `${context.subtitle}${context.scale !== 1 ? ' · значения ×10⁶' : ''}`;
      dom.matrixTable.setAttribute('aria-label', context.title);

      let html = '<thead><tr><th></th>' + context.labels.map((label, index) => `<th scope="col"${context.reductionPreview && fixed.includes(index) ? ' class="fixed-axis"' : ''}>${label}</th>`).join('') + '</tr></thead><tbody>';
      context.matrix.forEach((row, rowIndex) => {
        html += `<tr><th scope="row"${context.reductionPreview && fixed.includes(rowIndex) ? ' class="fixed-axis"' : ''}>${context.labels[rowIndex]}</th>`;
        row.forEach((value, columnIndex) => {
          const scaled = value / context.scale;
          const nonzero = Math.abs(value) > 1e-8;
          const selected = state.selectedCell && state.selectedCell.row === rowIndex && state.selectedCell.column === columnIndex;
          const currentCell = current.has(`${rowIndex}-${columnIndex}`);
          const eliminated = context.reductionPreview && (!free.includes(rowIndex) || !free.includes(columnIndex));
          if (eliminated) {
            html += `<td class="eliminated" aria-label="${context.labels[rowIndex]}, ${context.labels[columnIndex]}: исключено закреплением"><span class="eliminated-value">${formatNumber(scaled, compact ? 3 : 4)}</span></td>`;
          } else if (!nonzero) {
            html += '<td>0</td>';
          } else {
            const classes = ['matrix-cell'];
            if (context.kind === 'M') classes.push('mass');
            if (context.kind === 'Keff') classes.push('effective');
            if (currentCell) classes.push('current');
            if (selected) classes.push('selected');
            html += `<td><button class="${classes.join(' ')}" type="button" data-cell-row="${rowIndex}" data-cell-column="${columnIndex}" aria-label="${context.labels[rowIndex]}, ${context.labels[columnIndex]}: ${formatNumber(scaled, 5)}">${formatNumber(scaled, compact ? 3 : 4)}</button></td>`;
          }
        });
        html += '</tr>';
      });
      html += '</tbody>';
      dom.matrixTable.innerHTML = html;
    }

    function contributionAt(kind, globalRow, globalColumn, count = 2) {
      const contributions = [];
      for (let element = 0; element < count; element += 1) {
        const map = elements[element].map;
        const localRow = map.indexOf(globalRow);
        const localColumn = map.indexOf(globalColumn);
        if (localRow >= 0 && localColumn >= 0) {
          const value = elementMatrix(element, kind)[localRow][localColumn];
          if (Math.abs(value) > 1e-8) contributions.push({ element, value, localRow, localColumn });
        }
      }
      return contributions;
    }

    function renderCellDetail() {
      const context = displayContext();
      if (!state.selectedCell) {
        if (state.stage === 5) {
          dom.cellDetail.textContent = 'Нажмите K₄₄ или M₄₄: в общем узле 2 будет видна сумма вкладов e₁ и e₂.';
        } else if (state.stage === 6) {
          dom.cellDetail.textContent = 'Серые строки и столбцы принадлежат закреплённым СС и исключаются. Цветной блок [4, 5, 6] × [4, 5, 6] образует редуцированную матрицу 3 × 3.';
        } else if (state.stage >= 7) {
          dom.cellDetail.textContent = 'Нажмите ячейку редуцированной матрицы: будет показан её глобальный адрес и происхождение.';
        } else {
          dom.cellDetail.textContent = 'Выберите ненулевую ячейку, чтобы увидеть её происхождение.';
        }
        return;
      }

      const row = state.selectedCell.row;
      const column = state.selectedCell.column;
      const value = context.matrix[row][column];
      const symbol = context.kind === 'M' ? 'M' : context.kind === 'Keff' ? 'Keff' : 'K';

      if (context.kind === 'Keff') {
        const reduced = reducedMatrices();
        const a0 = 4 / state.timeStep ** 2;
        dom.cellDetail.textContent = `${symbol}[${row + 1},${column + 1}] = Kff + a₀Mff = ${formatNumber(reduced.K[row][column], 3)} + ${formatNumber(a0, 3)}·${formatNumber(reduced.M[row][column], 3)} = ${formatNumber(value, 3)}`;
        return;
      }

      if (state.stage === 1) {
        dom.cellDetail.textContent = `${symbol}local[${row + 1},${column + 1}] = ${formatPhysical(value, context.kind)}.`;
        return;
      }

      if (state.stage === 2 || state.stage === 3) {
        const element = state.stage - 2;
        dom.cellDetail.textContent = `${symbol}${element + 1}[${row + 1},${column + 1}] = (T${element + 1}ᵀ · ${symbol}local · T${element + 1})[${row + 1},${column + 1}] = ${formatPhysical(value, context.kind)}.`;
        return;
      }

      let globalRow = row;
      let globalColumn = column;
      if (state.stage >= 7) {
        globalRow = free[row];
        globalColumn = free[column];
      }
      const count = state.stage === 4 ? 1 : 2;
      const parts = contributionAt(context.kind, globalRow, globalColumn, count);
      const expression = parts.map(part => `${symbol}${part.element + 1}[${part.localRow + 1},${part.localColumn + 1}] = ${formatPhysical(part.value, context.kind)}`).join(' + ');
      dom.cellDetail.textContent = `${symbol}${globalRow + 1}${globalColumn + 1}: ${expression || 'нет ненулевых вкладов'}; итог = ${formatPhysical(value, context.kind)}.`;
    }

    function renderSummary() {
      const context = displayContext();
      const nonzeros = context.matrix.flat().filter(value => Math.abs(value) > 1e-8).length;
      let cards;
      if (state.stage === 6) {
        cards = [
          ['Закреплённые СС', '[1, 2, 3, 7, 8, 9]'],
          ['Свободные СС', '[4, 5, 6]'],
          ['Активный блок', '3 × 3 из 9 × 9']
        ];
      } else if (state.stage >= 7) {
        const reduced = reducedMatrices();
        cards = [
          ['Свободные СС', '[4, 5, 6]'],
          ['Kff[1,1]', `${formatNumber(reduced.K[0][0] / 1e6, 6)} ×10⁶`],
          ['Mff[1,1]', formatNumber(reduced.M[0][0], 3)]
        ];
      } else {
        cards = [
          ['Размер', `${context.matrix.length} × ${context.matrix.length}`],
          ['Ненулевых', `${nonzeros}`],
          ['Выбрано', state.matrix === 'M' ? 'матрица масс' : 'матрица жёсткости']
        ];
      }
      dom.summaryGrid.innerHTML = cards.map(([label, value]) => `<div class="summary-card"><span class="summary-card__label">${label}</span><span class="summary-card__value">${value}</span></div>`).join('');
    }

    function renderMatrixToggle() {
      document.querySelectorAll('[data-matrix]').forEach(button => {
        const isEffective = button.dataset.matrix === 'Keff';
        button.disabled = isEffective && state.stage < 7;
        button.setAttribute('aria-pressed', String(button.dataset.matrix === state.matrix));
      });
    }

    function renderNewmark() {
      const visible = state.stage === 7;
      dom.newmarkPanel.toggleAttribute('hidden', !visible);
      dom.timeStepValue.textContent = `${state.timeStep.toFixed(3)} s`;
      const a0 = 4 / state.timeStep ** 2;
      dom.effectiveFormula.textContent = `a₀ = 4/Δt² = ${formatNumber(a0, 1)};  Kₑff = Kff + a₀ · Mff`;
    }

    function render() {
      const stage = stages[state.stage];
      dom.stageTitle.textContent = stage.title;
      dom.stageDescription.textContent = stage.description;
      dom.stageBadge.textContent = `Шаг ${state.stage + 1} из ${stages.length}`;
      dom.equationChip.textContent = state.stage === 7 && state.matrix === 'Keff' ? 'Kₑff = Kff + 4/Δt² · Mff' : stage.equation;
      dom.prevButton.disabled = state.stage === 0;
      dom.nextButton.disabled = state.stage === stages.length - 1;
      dom.playButton.textContent = state.playing ? '❚❚ Пауза' : '▶ Запустить';

      renderStageStrip();
      renderElementPicker();
      renderModel();
      renderDofs();
      renderFormulaPanel();
      renderMatrixToggle();
      renderMatrix();
      renderCellDetail();
      renderSummary();
      renderNewmark();
    }

    function selectStage(index) {
      state.stage = Math.max(0, Math.min(stages.length - 1, index));
      state.selectedCell = null;
      if ([2, 4].includes(state.stage)) state.selectedElement = 0;
      if ([3, 5].includes(state.stage)) state.selectedElement = 1;
      if (state.stage < 7 && state.matrix === 'Keff') state.matrix = 'K';
      if (state.stage === 7) state.matrix = 'Keff';
      if (state.stage === stages.length - 1 && state.playing) stopPlayback();
      render();
    }

    function selectElement(index) {
      state.selectedElement = index;
      if (state.stage === 0) state.stage = 1;
      state.selectedCell = null;
      render();
    }

    function selectMatrix(kind) {
      if (kind === 'Keff' && state.stage < 7) return;
      state.matrix = kind;
      state.selectedCell = null;
      render();
    }

    function stopPlayback() {
      state.playing = false;
      window.clearInterval(state.timer);
      state.timer = null;
    }

    function restartPlaybackTimer() {
      window.clearInterval(state.timer);
      const delays = { 1: 2800, 2: 1800, 3: 1100 };
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

    dom.prevButton.addEventListener('click', () => selectStage(state.stage - 1));
    dom.nextButton.addEventListener('click', () => selectStage(state.stage + 1));
    dom.playButton.addEventListener('click', togglePlayback);
    dom.resetButton.addEventListener('click', () => {
      stopPlayback();
      state.stage = 0;
      state.matrix = 'K';
      state.selectedElement = 0;
      state.selectedCell = null;
      render();
    });

    dom.speedRange.addEventListener('input', () => { if (state.playing) restartPlaybackTimer(); });

    dom.timeStepRange.addEventListener('input', event => {
      state.timeStep = Number(event.target.value);
      state.selectedCell = null;
      render();
    });

    document.addEventListener('click', event => {
      const stageButton = event.target.closest('[data-stage]');
      if (stageButton) {
        stopPlayback();
        selectStage(Number(stageButton.dataset.stage));
        return;
      }
      const matrixButton = event.target.closest('[data-matrix]');
      if (matrixButton) {
        selectMatrix(matrixButton.dataset.matrix);
        return;
      }
      const elementButton = event.target.closest('[data-element-button]');
      if (elementButton) {
        selectElement(Number(elementButton.dataset.elementButton));
        return;
      }
      const elementGroup = event.target.closest('[data-element]');
      if (elementGroup) {
        selectElement(Number(elementGroup.dataset.element));
        return;
      }
      const cell = event.target.closest('[data-cell-row]');
      if (cell) {
        state.selectedCell = { row: Number(cell.dataset.cellRow), column: Number(cell.dataset.cellColumn) };
        renderMatrix();
        renderCellDetail();
      }
    });

    document.querySelectorAll('[data-element]').forEach(group => {
      group.addEventListener('keydown', event => {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          selectElement(Number(group.dataset.element));
        }
      });
    });

    render();
