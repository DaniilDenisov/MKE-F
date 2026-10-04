'use strict';

    const stages = [
      { short: 'Контейнеры', title: 'Два контейнера для одной системы', description: 'Плотная матрица заранее выделяет n² ячеек. Триплетная сборка резервирует линейные массивы только под вклады элементов.' },
      { short: 'Один элемент', title: 'Карта степеней свободы создаёт адреса', description: 'Выберите элемент: ndgrid создаёт четыре пары адресов, а их цвета сразу связывают строки, столбцы и коэффициенты полной матрицы.' },
      { short: 'Триплеты', title: 'Все элементы записаны последовательно', description: 'Выберите элемент: четыре коэффициента и номера задействованных строк и столбцов полной матрицы получат цвета соответствующих триплетов.' },
      { short: 'sparse()', title: 'GNU Octave объединяет повторные адреса', description: 'sparse(rows, columns, values, n, n) суммирует совпадающие координаты и создаёт глобальную разреженную матрицу.' },
      { short: 'Решение', title: 'Разреженная матрица поступает в решатель', description: 'После исключения закреплённых СС оператор обратной косой черты решает систему; вручную превращать Kff или Mff в полную матрицу не требуется.' }
    ];

    const state = { stage: 0, elementCount: 8, matrixType: 'K', selectedElement: 0 };

    const dom = {
      stageStrip: document.getElementById('stageStrip'),
      stageTitle: document.getElementById('stageTitle'),
      stageDescription: document.getElementById('stageDescription'),
      stageBadge: document.getElementById('stageBadge'),
      prevButton: document.getElementById('prevButton'),
      nextButton: document.getElementById('nextButton'),
      resetButton: document.getElementById('resetButton'),
      elementRange: document.getElementById('elementRange'),
      elementCountValue: document.getElementById('elementCountValue'),
      chain: document.getElementById('chain'),
      matrixTable: document.getElementById('matrixTable'),
      denseChip: document.getElementById('denseChip'),
      sparseChip: document.getElementById('sparseChip'),
      denseCode: document.getElementById('denseCode'),
      dofsArray: document.getElementById('dofsArray'),
      elementRowsArray: document.getElementById('elementRowsArray'),
      elementColumnsArray: document.getElementById('elementColumnsArray'),
      elementValuesName: document.getElementById('elementValuesName'),
      elementValuesArray: document.getElementById('elementValuesArray'),
      elementValuesNote: document.getElementById('elementValuesNote'),
      tripletList: document.getElementById('tripletList'),
      sparseCall: document.getElementById('sparseCall'),
      sparseExplanation: document.getElementById('sparseExplanation'),
      uniqueCount: document.getElementById('uniqueCount'),
      duplicateCount: document.getElementById('duplicateCount'),
      fillRatio: document.getElementById('fillRatio'),
      productionCode: document.getElementById('productionCode'),
      denseMetric: document.getElementById('denseMetric'),
      bufferMetric: document.getElementById('bufferMetric'),
      sparseMetric: document.getElementById('sparseMetric'),
      denseBar: document.getElementById('denseBar'),
      bufferBar: document.getElementById('bufferBar'),
      sparseBar: document.getElementById('sparseBar'),
      solveTitle: document.getElementById('solveTitle'),
      solveText: document.getElementById('solveText'),
      solveEquation: document.getElementById('solveEquation')
    };

    function formatNumber(value) {
      if (Math.abs(value) < 1e-10) return '0';
      return String(Number(value.toFixed(3))).replace('-', '−');
    }

    function elementData(index) {
      const dofs = [index, index + 1];
      if (state.matrixType === 'K') {
        const stiffness = 10 + 5 * (index % 3);
        return { dofs, values: [stiffness, -stiffness, -stiffness, stiffness], symbol: `k${index + 1}` };
      }
      const mass = 6 + 3 * (index % 3);
      return { dofs, values: [mass / 3, mass / 6, mass / 6, mass / 3], symbol: `m${index + 1}` };
    }

    function elementTriplets(index) {
      const data = elementData(index);
      return [
        { row: data.dofs[0], column: data.dofs[0], value: data.values[0], element: index, slot: 0 },
        { row: data.dofs[1], column: data.dofs[0], value: data.values[1], element: index, slot: 1 },
        { row: data.dofs[0], column: data.dofs[1], value: data.values[2], element: index, slot: 2 },
        { row: data.dofs[1], column: data.dofs[1], value: data.values[3], element: index, slot: 3 }
      ];
    }

    function visibleElements() {
      if (state.stage === 0) return [];
      if (state.stage === 1) return [state.selectedElement];
      return Array.from({ length: state.elementCount }, (_, index) => index);
    }

    function rawTriplets() {
      return visibleElements().flatMap(elementTriplets);
    }

    function allTriplets() {
      return Array.from({ length: state.elementCount }, (_, index) => index).flatMap(elementTriplets);
    }

    function coalesce(triplets) {
      const entries = new Map();
      triplets.forEach(item => {
        const key = `${item.row}-${item.column}`;
        entries.set(key, (entries.get(key) || 0) + item.value);
      });
      return entries;
    }

    function matrixFrom(triplets) {
      const size = state.elementCount + 1;
      const matrix = Array.from({ length: size }, () => Array(size).fill(0));
      triplets.forEach(item => { matrix[item.row][item.column] += item.value; });
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
        if (Math.abs(divisor) < 1e-12) return [];
        for (let column = pivot; column <= n; column += 1) augmented[pivot][column] /= divisor;
        for (let row = 0; row < n; row += 1) {
          if (row === pivot) continue;
          const factor = augmented[row][pivot];
          for (let column = pivot; column <= n; column += 1) augmented[row][column] -= factor * augmented[pivot][column];
        }
      }
      return augmented.map(row => row[n]);
    }

    function renderStages() {
      dom.stageStrip.innerHTML = stages.map((stage, index) => `
        <button class="stage-button${index < state.stage ? ' is-complete' : ''}" type="button" data-stage="${index}" ${index === state.stage ? 'aria-current="step"' : ''}>
          <span>${index < state.stage ? '✓ завершено' : `шаг ${index + 1}`}</span>${stage.short}
        </button>
      `).join('');
    }

    function renderChain() {
      let html = '<span class="node is-fixed" data-label="u₁=0" aria-hidden="true"></span>';
      for (let index = 0; index < state.elementCount; index += 1) {
        html += `<button class="member" type="button" data-element="${index}" data-label="e${index + 1}" aria-label="Выбрать элемент ${index + 1}" aria-pressed="${index === state.selectedElement}"></button>`;
        html += `<span class="node" data-label="u${index + 2}" aria-hidden="true"></span>`;
      }
      dom.chain.innerHTML = html;
    }

    function renderMatrix() {
      const size = state.elementCount + 1;
      const matrix = matrixFrom(rawTriplets());
      const selectedTriplets = elementTriplets(state.selectedElement);
      const currentCells = new Map(selectedTriplets.map(item => [`${item.row}-${item.column}`, item]));
      const selectedDOFs = elementData(state.selectedElement).dofs;
      const showSelectedColors = state.stage === 1 || state.stage === 2;
      let html = '<thead><tr><th></th>';
      for (let column = 0; column < size; column += 1) {
        let label = String(column + 1);
        if (showSelectedColors && column === selectedDOFs[0]) label = `<span class="axis-number axis-pair-0-1">${label}</span>`;
        if (showSelectedColors && column === selectedDOFs[1]) label = `<span class="axis-number axis-pair-2-3">${label}</span>`;
        html += `<th scope="col">${label}</th>`;
      }
      html += '</tr></thead><tbody>';
      matrix.forEach((row, rowIndex) => {
        let rowLabel = String(rowIndex + 1);
        if (showSelectedColors && rowIndex === selectedDOFs[0]) rowLabel = `<span class="axis-number axis-pair-0-2">${rowLabel}</span>`;
        if (showSelectedColors && rowIndex === selectedDOFs[1]) rowLabel = `<span class="axis-number axis-pair-1-3">${rowLabel}</span>`;
        html += `<tr><th scope="row">${rowLabel}</th>`;
        row.forEach((value, columnIndex) => {
          const classes = [Math.abs(value) > 1e-10 ? 'nonzero' : 'zero'];
          const selectedTriplet = currentCells.get(`${rowIndex}-${columnIndex}`);
          if (showSelectedColors && selectedTriplet) classes.push('current', `slot-${selectedTriplet.slot}`);
          const contribution = showSelectedColors && selectedTriplet
            ? `; вклад выбранного элемента = ${formatNumber(selectedTriplet.value)}`
            : '';
          html += `<td class="${classes.join(' ')}" title="(${rowIndex + 1}, ${columnIndex + 1}) = ${formatNumber(value)}${contribution}">${size > 12 ? '•' : formatNumber(value)}</td>`;
        });
        html += '</tr>';
      });
      html += '</tbody>';
      dom.matrixTable.classList.toggle('compact', size > 12);
      dom.matrixTable.innerHTML = html;
      dom.matrixTable.setAttribute('aria-label', `Плотное представление матрицы ${state.matrixType} размером ${size} на ${size}`);
      dom.denseChip.textContent = `${size * size} ячеек`;
    }

    function renderElementArrays() {
      const data = elementData(state.selectedElement);
      const rows = [data.dofs[0], data.dofs[1], data.dofs[0], data.dofs[1]].map(value => value + 1);
      const columns = [data.dofs[0], data.dofs[0], data.dofs[1], data.dofs[1]].map(value => value + 1);
      const coloredValues = values => `[${values.map((value, slot) => `<span class="value-token slot-${slot}">${formatNumber(value)}</span>`).join(', ')}]`;
      dom.dofsArray.innerHTML = `[${data.dofs.map(value => `<span class="value-token value-token--dof">${value + 1}</span>`).join(', ')}]`;
      dom.elementRowsArray.innerHTML = coloredValues(rows);
      dom.elementColumnsArray.innerHTML = coloredValues(columns);
      dom.elementValuesName.textContent = state.matrixType === 'K' ? 'stiffness(:)' : 'mass(:)';
      dom.elementValuesArray.innerHTML = coloredValues(data.values);
      dom.elementValuesNote.textContent = state.matrixType === 'K' ? 'коэффициенты жёсткости' : 'коэффициенты массы';
    }

    function renderTriplets() {
      const triplets = rawTriplets();
      const addressCounts = new Map();
      triplets.forEach(item => {
        const key = `${item.row}-${item.column}`;
        addressCounts.set(key, (addressCounts.get(key) || 0) + 1);
      });
      const visible = triplets.slice(0, 56);
      dom.tripletList.innerHTML = visible.length ? visible.map(item => {
        const key = `${item.row}-${item.column}`;
        const classes = ['triplet', `slot-${item.slot}`];
        if (item.element === state.selectedElement) classes.push('is-current');
        if ((addressCounts.get(key) || 0) > 1) classes.push('is-duplicate');
        return `<span class="${classes.join(' ')}">(${item.row + 1}, ${item.column + 1}, ${formatNumber(item.value)})</span>`;
      }).join('') + (triplets.length > visible.length ? `<span class="triplet">… ещё ${triplets.length - visible.length}</span>` : '') : '<span class="triplet">буферы пока пусты</span>';
      dom.sparseChip.textContent = `${triplets.length} вкладов-триплетов`;
    }

    function renderSparseSummary() {
      const size = state.elementCount + 1;
      const triplets = state.stage < 2 ? rawTriplets() : allTriplets();
      const unique = coalesce(triplets);
      const duplicates = triplets.length - unique.size;
      dom.uniqueCount.textContent = String(unique.size);
      dom.duplicateCount.textContent = String(duplicates);
      dom.fillRatio.textContent = `${((unique.size / (size * size)) * 100).toFixed(1)} %`;
      const target = state.matrixType === 'K' ? 'globalK' : 'globalM';
      const values = state.matrixType === 'K' ? 'stiffnessValues' : 'massValues';
      dom.sparseCall.textContent = `${target} = sparse(rows, columns, ${values}, n, n);`;
      dom.sparseExplanation.textContent = state.stage >= 3
        ? `GNU Octave объединил ${duplicates} повторяющихся адресов. Нулевые позиции полной матрицы не нужно перечислять при последующих операциях.`
        : 'Библиотечная функция GNU Octave принимает координаты и значения, складывает записи с одинаковым адресом и создаёт разреженную матрицу заданного размера.';
    }

    function renderCode() {
      const size = state.elementCount + 1;
      const matrixName = state.matrixType;
      const lower = state.matrixType === 'K' ? 'stiffness' : 'mass';
      const denseSnippets = [
        `<span class="code-comment">% Плотный квадрат существует целиком</span>\n${matrixName} = zeros(${size}, ${size});`,
        `dofs = elements(${state.selectedElement + 1}).dofs;\n${matrixName}(dofs,dofs) = ${matrixName}(dofs,dofs) + element.${lower};`,
        `for e = 1:${state.elementCount}\n  dofs = elements(e).dofs;\n  ${matrixName}(dofs,dofs) += elements(e).${lower};\nend`,
        `<span class="code-comment">% Результат тот же, но нули занимают ячейки</span>\n${matrixName} = full(global${matrixName});`,
        `freeDOFs = 2:${size};\n${matrixName}ff = ${matrixName}(freeDOFs,freeDOFs);`
      ];
      const productionSnippets = [
        `entryCount = Σ numel(dofs)²;\nrows = zeros(entryCount,1);\ncolumns = zeros(entryCount,1);`,
        `[elementRows,elementColumns] = ndgrid(dofs,dofs);\nrows(indices) = elementRows(:);\ncolumns(indices) = elementColumns(:);`,
        `${lower}Values(indices) = element.${lower}(:);\ncursor = cursor + numel(dofs)^2;`,
        `<span class="code-accent">global${matrixName} = sparse(rows,columns,${lower}Values,${size},${size});</span>`,
        state.matrixType === 'K'
          ? `reducedK = model.stiffness(freeDOFs,freeDOFs);\nu(freeDOFs) = reducedK \\ F(freeDOFs);`
          : `reducedM = model.mass(freeDOFs,freeDOFs);\na0 = reducedM \\ (F0 - reducedK*u0);`
      ];
      dom.denseCode.innerHTML = denseSnippets[state.stage];
      dom.productionCode.innerHTML = productionSnippets[state.stage];
    }

    function renderMetrics() {
      const size = state.elementCount + 1;
      const denseWords = 2 * size * size;
      const bufferWords = 16 * state.elementCount;
      const nonzeroWords = 2 * coalesce(allTriplets()).size;
      const maximum = Math.max(denseWords, bufferWords, nonzeroWords, 1);
      dom.denseMetric.textContent = `${denseWords} чисел`;
      dom.bufferMetric.textContent = `${bufferWords} чисел`;
      dom.sparseMetric.textContent = `${nonzeroWords} значений`;
      dom.denseBar.style.width = `${100 * denseWords / maximum}%`;
      dom.bufferBar.style.width = `${100 * bufferWords / maximum}%`;
      dom.sparseBar.style.width = `${100 * nonzeroWords / maximum}%`;
    }

    function renderSolve() {
      const matrix = matrixFrom(allTriplets());
      const reduced = matrix.slice(1).map(row => row.slice(1));
      const rightHandSide = Array(state.elementCount).fill(0);
      rightHandSide[rightHandSide.length - 1] = 30;
      const solution = solveLinearSystem(reduced, rightHandSide);
      const last = solution.length ? formatNumber(solution[solution.length - 1]) : '—';
      if (state.matrixType === 'K') {
        dom.solveTitle.textContent = 'Статическая система';
        dom.solveText.textContent = `Первая СС закреплена, Fₙ = 30. Для текущих коэффициентов получено uₙ = ${last}.`;
        dom.solveEquation.textContent = 'u(freeDOFs) = reducedK \\ loads(freeDOFs)';
      } else {
        dom.solveTitle.textContent = 'Начальное ускорение';
        dom.solveText.textContent = `Та же операция применяется к матрице масс: при условном правом столбце 30 получено aₙ = ${last}.`;
        dom.solveEquation.textContent = 'a₀ = reducedM \\ (F₀ − reducedK·u₀)';
      }
    }

    function render() {
      const stage = stages[state.stage];
      dom.stageTitle.textContent = stage.title;
      dom.stageDescription.textContent = stage.description;
      dom.stageBadge.textContent = `Шаг ${state.stage + 1} из ${stages.length}`;
      dom.prevButton.disabled = state.stage === 0;
      dom.nextButton.disabled = state.stage === stages.length - 1;
      dom.elementCountValue.textContent = String(state.elementCount);
      document.querySelectorAll('[data-matrix]').forEach(button => button.setAttribute('aria-pressed', String(button.dataset.matrix === state.matrixType)));
      renderStages();
      renderChain();
      renderMatrix();
      renderElementArrays();
      renderTriplets();
      renderSparseSummary();
      renderCode();
      renderMetrics();
      renderSolve();
    }

    dom.prevButton.addEventListener('click', () => { if (state.stage > 0) { state.stage -= 1; render(); } });
    dom.nextButton.addEventListener('click', () => { if (state.stage < stages.length - 1) { state.stage += 1; render(); } });
    dom.resetButton.addEventListener('click', () => { state.stage = 0; state.selectedElement = 0; render(); });
    dom.elementRange.addEventListener('input', event => {
      state.elementCount = Number(event.target.value);
      state.selectedElement = Math.min(state.selectedElement, state.elementCount - 1);
      render();
    });
    document.addEventListener('click', event => {
      const stageButton = event.target.closest('[data-stage]');
      if (stageButton) { state.stage = Number(stageButton.dataset.stage); render(); return; }
      const matrixButton = event.target.closest('[data-matrix]');
      if (matrixButton) { state.matrixType = matrixButton.dataset.matrix; render(); return; }
      const elementButton = event.target.closest('[data-element]');
      if (elementButton) { state.selectedElement = Number(elementButton.dataset.element); if (state.stage === 0) state.stage = 1; render(); }
    });

    render();
