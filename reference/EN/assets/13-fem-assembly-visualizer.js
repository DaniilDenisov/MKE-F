'use strict';

    const stages = [
      {
        short: 'Model',
        title: 'Model and degrees of freedom',
        description: 'Each node has one axial DOF. Numbering determines coefficient addresses in the global matrix.',
        equation: 'K = 0'
      },
      {
        short: 'Element',
        title: 'Element matrix',
        description: 'Select any member: its matrix depends on k = EA/L, and its local numbers map to two global DOFs.',
        equation: 'k⁽ᵉ⁾ = (EA/L) · [1 −1; −1 1]'
      },
      {
        short: 'Contribution e₁',
        title: 'Assembling element e₁',
        description: 'The map [u₁, u₂] places the first element’s four coefficients in K.',
        equation: 'K([1,2],[1,2]) += k⁽¹⁾'
      },
      {
        short: 'Contribution e₂',
        title: 'Assembling element e₂',
        description: 'The second matrix enters rows and columns [u₂, u₃]. Two contributions add at K₂₂.',
        equation: 'K([2,3],[2,3]) += k⁽²⁾'
      },
      {
        short: 'Contribution e₃',
        title: 'Assembling element e₃',
        description: 'After the third contribution, the global matrix is fully assembled with a tridiagonal pattern.',
        equation: 'K([3,4],[3,4]) += k⁽³⁾'
      },
      {
        short: 'Fixed support',
        title: 'Fixed support u₁ = 0',
        description: 'The constrained DOF is removed. The free block Kff is solved for [u₂, u₃, u₄].',
        equation: 'Kff · uf = Ff'
      },
      {
        short: 'Solution',
        title: 'Displacements and reaction',
        description: 'Free displacements are expanded to the full vector; the reaction comes from the original system: R = Ku − F.',
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
        if (Math.abs(divisor) < 1e-12) throw new Error('Singular matrix');
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
          <span>${index < state.stage ? '✓ complete' : `step ${index + 1}`}</span>
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
      dom.localTitle.textContent = `Element e${element + 1}`;
      dom.localStiffness.textContent = `k${element + 1} = E${element + 1}A${element + 1}/L${element + 1} = ${formatNumber(k)} kN/mm`;
      dom.dofMapChip.textContent = `global [u${element + 1}, u${element + 2}]`;
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

      dom.loadLabel.textContent = `P = ${formatNumber(state.load)} kN`;
      const solved = state.stage === 6;
      dom.deformedShape.toggleAttribute('hidden', !solved);
      dom.reactionArrow.toggleAttribute('hidden', !solved);
      dom.reactionLabel.textContent = `R₁ = −${formatNumber(state.load)} kN`;

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
      dom.globalMatrix.setAttribute('aria-label', state.stage >= 5 ? 'Reduced stiffness matrix' : 'Global stiffness matrix');
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
        : '<span class="triplet">empty so far</span>';

      const matrix = displayedMatrix();
      const csr = matrixToCsr(matrix);
      dom.csrTitle.textContent = `After sparse() · ${matrix.length} rows · educational row-wise view`;
      dom.csrValues.textContent = `[${csr.values.map(value => formatNumber(value)).join(', ')}]`;
      dom.csrColumns.textContent = `[${csr.columns.join(', ')}]`;
      dom.csrRowPtr.textContent = `[${csr.rowPtr.join(', ')}]`;
      const denseNumbers = matrix.length * matrix.length;
      dom.storageComparison.textContent = `Dense: ${denseNumbers} numbers · sparse: ${csr.values.length} nonzero values; columns and rowPtr below are an educational row-wise representation only`;
    }

    function operationText() {
      const k = state.stiffness;
      const P = state.load;
      if (state.stage === 0) return 'K is currently zero. The element map determines which four cells receive its contribution.';
      if (state.stage === 1) return `For selected e${state.selectedElement + 1}: k = EA/L = ${formatNumber(k[state.selectedElement])} kN/mm; local coefficients have not yet been added to K.`;
      if (state.stage === 2) return `K₁₁ += k₁; K₁₂ += −k₁; K₂₁ += −k₁; K₂₂ += k₁. Result: K₂₂ = ${formatNumber(k[0])}.`;
      if (state.stage === 3) return `Repeated address is summed: K₂₂ = k₁ + k₂ = ${formatNumber(k[0])} + ${formatNumber(k[1])} = ${formatNumber(k[0] + k[1])}.`;
      if (state.stage === 4) return `K₃₃ = k₂ + k₃ = ${formatNumber(k[1])} + ${formatNumber(k[2])} = ${formatNumber(k[1] + k[2])}. Global K is assembled.`;
      if (state.stage === 5) return `u₁ = 0, so row and column 1 are excluded from Kff. Ff = [0, 0, ${formatNumber(P)}]ᵀ kN; original K is retained for reactions.`;
      const u = solution().u;
      return `u = [${u.map(value => formatNumber(value)).join(', ')}]ᵀ mm; R₁ = −${formatNumber(P)} kN; each element carries N = ${formatNumber(P)} kN.`;
    }

    function renderResults() {
      const count = assembledElementCount();
      if (state.stage === 6) {
        const result = solution();
        dom.resultGrid.innerHTML = `
          <div class="result-card">
            <span class="result-card__label">Displacements u, mm</span>
            <span class="result-card__value">[${result.u.map(value => formatNumber(value)).join(', ')}]ᵀ</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">Reaction R₁</span>
            <span class="result-card__value">−${formatNumber(state.load)} kN</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">Axial forces</span>
            <span class="result-card__value">N₁ = N₂ = N₃ = ${formatNumber(state.load)} kN</span>
          </div>
        `;
      } else {
        dom.resultGrid.innerHTML = `
          <div class="result-card">
            <span class="result-card__label">Elements assembled</span>
            <span class="result-card__value">${count} / 3</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">System size to solve</span>
            <span class="result-card__value">${state.stage >= 5 ? '3 × 3 (from 4 × 4)' : '4 × 4'}</span>
          </div>
          <div class="result-card">
            <span class="result-card__label">State</span>
            <span class="result-card__value">${state.stage === 0 ? 'K not assembled' : state.stage === 1 ? 'Element matrix' : state.stage < 4 ? 'Assembling' : state.stage === 4 ? 'K assembled' : 'Kff ready'}</span>
          </div>
        `;
      }

      dom.footnote.textContent = state.stage >= 5
        ? 'Reduction removes the constrained DOF. Full K is needed to recover the reaction R = Ku − F.'
        : 'The dense view shows every cell, including zeros. The program passes triplets to GNU Octave sparse(), which sums repeated addresses and creates a sparse matrix.';
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
      dom.stageBadge.textContent = `Step ${state.stage + 1} of ${stages.length}`;
      dom.equationChip.textContent = stage.equation;
      dom.matrixTitle.textContent = state.stage >= 5 ? 'Reduced matrix Kff' : 'Global stiffness matrix K';
      dom.operationLine.textContent = operationText();
      dom.prevButton.disabled = state.stage === 0;
      dom.nextButton.disabled = state.stage === stages.length - 1;
      dom.playButton.textContent = state.playing ? '❚❚ Pause' : '▶ Play';

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
        document.getElementById(`stiffnessValue${index + 1}`).textContent = `${formatNumber(state.stiffness[index])} kN/mm`;
        render();
      });
    });

    document.getElementById('loadRange').addEventListener('input', event => {
      state.load = Number(event.target.value);
      document.getElementById('loadValue').textContent = `${formatNumber(state.load)} kN`;
      render();
    });

    render();
