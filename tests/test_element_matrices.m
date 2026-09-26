function test_element_matrices()
%TEST_ELEMENT_MATRICES Проверка матриц плоских стержневых и рамных элементов.
% Эталонные матрицы независимо вычисляются ниже и сравниваются с матрицами,
% полученными через открытый интерфейс сборки каждого элемента.

testTrussOrientations();
testBeamOrientations();
testInvalidElementData();
end

function testTrussOrientations()
area = 0.02;
youngsModulus = 210e9;
density = 7850;
properties = [area youngsModulus density];

% Все три элемента имеют длину 5. Одинаковые длина и свойства позволяют напрямую
% сравнивать спектры матриц после поворота элемента.
coordinates = {
    [0 0; 5 0]
    [0 0; 0 5]
    [1 2; 4 6]
};
referenceStiffnessSpectrum = [];
referenceMassSpectrum = [];

for i = 1:numel(coordinates)
    coords = coordinates{i};
    [K, M] = assembleTruss(coords, properties);
    [expectedK, expectedM, length] = expectedTrussMatrices(coords, properties);

    assertClose(K, expectedK, 1e-12, 1e-8, ...
        'Truss stiffness matrix has incorrect coefficients.');
    assertClose(M, expectedM, 1e-12, 1e-12, ...
        'Truss mass matrix has incorrect coefficients.');
    assertSymmetric(K, 'Truss stiffness matrix is not symmetric.');
    assertSymmetric(M, 'Truss mass matrix is not symmetric.');

    % У плоского двухузлового стержня одна деформационная форма. Остальные три
    % формы — движения с нулевой энергией, поэтому ранг матрицы жёсткости равен 1.
    assertNumericalRank(K, 1, ...
        'A free truss element must have three zero-energy modes.');

    % При единичном переносе квадратичная форма должна давать физическую массу
    % элемента rho*A*L для каждого глобального направления.
    rigidX = [1; 0; 1; 0];
    rigidY = [0; 1; 0; 1];
    expectedMass = density * area * length;
    assertClose(rigidX.' * M * rigidX, expectedMass, 1e-12, 1e-12, ...
        'Truss total mass in global X is incorrect.');
    assertClose(rigidY.' * M * rigidY, expectedMass, 1e-12, 1e-12, ...
        'Truss total mass in global Y is incorrect.');

    % Ортогональный поворот координат сохраняет собственные значения матрицы.
    % Предварительная симметризация исключает влияние ошибок округления.
    stiffnessSpectrum = sort(eig((K + K.') / 2));
    massSpectrum = sort(eig((M + M.') / 2));
    if isempty(referenceStiffnessSpectrum)
        referenceStiffnessSpectrum = stiffnessSpectrum;
        referenceMassSpectrum = massSpectrum;
    else
        assertClose(stiffnessSpectrum, referenceStiffnessSpectrum, ...
            1e-12, 1e-6, 'Truss stiffness is not rotation invariant.');
        assertClose(massSpectrum, referenceMassSpectrum, ...
            1e-12, 1e-12, 'Truss mass is not rotation invariant.');
    end
end
end

function testBeamOrientations()
area = 0.02;
youngsModulus = 210e9;
density = 7850;
momentOfInertia = 8e-6;
properties = [area youngsModulus density momentOfInertia];

% Горизонтальный, вертикальный и наклонный элемент 3-4-5 проверяют все члены
% преобразования при одинаковой длине и одинаковых физических свойствах.
coordinates = {
    [0 0; 5 0]
    [0 0; 0 5]
    [1 2; 4 6]
};
referenceStiffnessSpectrum = [];
referenceMassSpectrum = [];

for i = 1:numel(coordinates)
    coords = coordinates{i};
    [K, M] = assembleBeam(coords, properties);
    [expectedK, expectedM, length] = expectedBeamMatrices(coords, properties);

    assertClose(K, expectedK, 1e-12, 1e-8, ...
        'Beam stiffness matrix has incorrect coefficients.');
    assertClose(M, expectedM, 1e-12, 1e-12, ...
        'Beam mass matrix has incorrect coefficients.');
    assertSymmetric(K, 'Beam stiffness matrix is not symmetric.');
    assertSymmetric(M, 'Beam mass matrix is not symmetric.');

    % Свободный плоский рамный элемент имеет три деформационные формы и три формы
    % движения как твёрдого тела, поэтому ранг матрицы жёсткости 6x6 равен 3.
    assertNumericalRank(K, 3, ...
        'A free beam element must have three rigid-body modes.');

    % При поступательном движении вращательные степени свободы равны нулю.
    % Поэтому квадратичная форма масс должна давать rho*A*L по обеим глобальным осям.
    rigidX = [1; 0; 0; 1; 0; 0];
    rigidY = [0; 1; 0; 0; 1; 0];
    expectedMass = density * area * length;
    assertClose(rigidX.' * M * rigidX, expectedMass, 1e-12, 1e-12, ...
        'Beam total mass in global X is incorrect.');
    assertClose(rigidY.' * M * rigidY, expectedMass, 1e-12, 1e-12, ...
        'Beam total mass in global Y is incorrect.');

    % Одинаковые спектры подтверждают, что изменение только ориентации элемента
    % не изменяет его физическую жёсткость и инерцию.
    stiffnessSpectrum = sort(eig((K + K.') / 2));
    massSpectrum = sort(eig((M + M.') / 2));
    if isempty(referenceStiffnessSpectrum)
        referenceStiffnessSpectrum = stiffnessSpectrum;
        referenceMassSpectrum = massSpectrum;
    else
        assertClose(stiffnessSpectrum, referenceStiffnessSpectrum, ...
            1e-12, 1e-5, 'Beam stiffness is not rotation invariant.');
        assertClose(massSpectrum, referenceMassSpectrum, ...
            1e-12, 1e-10, 'Beam mass is not rotation invariant.');
    end
end
end

function testInvalidElementData()
validCoords = [0 0; 1 0];

% Вырожденные и неконечные координаты должны быть отклонены до того, как
% преобразование выполнит деление на длину элемента.
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupTruss([0 0; 0 0], [1 2e11 7850]));
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupBeam([0 0; 0 0], [1 2e11 7850 1]));
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupTruss([0 0; Inf 0], [1 2e11 7850]));
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupBeam([0 0; NaN 1], [1 2e11 7850 1]));

% Неположительные или неконечные физические свойства лишают модель смысла
% и могут скрыть дальнейшие вырождения, поэтому отклоняем их при настройке.
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupTruss(validCoords, [0 2e11 7850]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupTruss(validCoords, [1 -2e11 7850]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupTruss(validCoords, [1 2e11 0]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupBeam(validCoords, [1 2e11 7850 0]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupBeam(validCoords, [1 2e11 Inf 1]));
end

function [K, M] = assembleTruss(coords, properties)
% Единый структурный интерфейс предоставляет матрицы элемента напрямую.
element = createStructuralElement(112, coords, [1 2], properties, ...
    [1 2; 3 4]);
K = element.stiffness;
M = element.mass;
end

function [K, M] = assembleBeam(coords, properties)
element = createStructuralElement(113, coords, [1 2], properties, ...
    [1 2 3; 4 5 6]);
K = element.stiffness;
M = element.mass;
end

function setupTruss(coords, properties)
createStructuralElement(112, coords, [1 2], properties, [1 2; 3 4]);
end

function setupBeam(coords, properties)
createStructuralElement(113, coords, [1 2], properties, ...
    [1 2 3; 4 5 6]);
end

function [K, M, length] = expectedTrussMatrices(coords, properties)
% Стандартные глобальные матрицы жёсткости и согласованных масс плоского стержня.
delta = coords(2,:) - coords(1,:);
length = norm(delta);
c = delta(1) / length;
s = delta(2) / length;
area = properties(1);
youngsModulus = properties(2);
density = properties(3);

K = area * youngsModulus / length * ...
    [ c*c  c*s -c*c -c*s; ...
      c*s  s*s -c*s -s*s; ...
     -c*c -c*s  c*c  c*s; ...
     -c*s -s*s  c*s  s*s];
M = density * area * length / 6 * ...
    [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
end

function [K, M, length] = expectedBeamMatrices(coords, properties)
% Стандартные локальные матрицы рамного элемента Эйлера-Бернулли независимо
% преобразуются в глобальные координаты для выявления ошибок знаков и коэффициентов.
delta = coords(2,:) - coords(1,:);
length = norm(delta);
c = delta(1) / length;
s = delta(2) / length;
area = properties(1);
youngsModulus = properties(2);
density = properties(3);
momentOfInertia = properties(4);

axial = youngsModulus * area / length;
bending12 = 12 * youngsModulus * momentOfInertia / length^3;
bending6 = 6 * youngsModulus * momentOfInertia / length^2;
bending4 = 4 * youngsModulus * momentOfInertia / length;
bending2 = 2 * youngsModulus * momentOfInertia / length;
localK = [ axial       0          0 -axial        0          0; ...
               0 bending12   bending6      0 -bending12   bending6; ...
               0  bending6   bending4      0  -bending6   bending2; ...
          -axial       0          0  axial        0          0; ...
               0 -bending12 -bending6      0  bending12  -bending6; ...
               0  bending6   bending2      0  -bending6   bending4];
localM = density * area * length / 420 * ...
    [140       0            0  70       0            0; ...
       0     156    22*length   0      54   -13*length; ...
       0 22*length 4*length^2   0 13*length -3*length^2; ...
      70       0            0 140       0            0; ...
       0      54    13*length   0     156   -22*length; ...
       0 -13*length -3*length^2 0 -22*length 4*length^2];
transform = [ c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; ...
              0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
K = transform.' * localK * transform;
M = transform.' * localM * transform;
end

function assertSymmetric(matrix, message)
% Сравниваем нарушение симметрии с масштабом матрицы, а не только с абсолютным допуском.
errorNorm = norm(matrix - matrix.', inf);
scale = max(norm(matrix, inf), 1);
if errorNorm > 1e-12 * scale
    error('MKEF:VerificationFailed', '%s Error=%g, scale=%g.', ...
        message, errorNorm, scale);
end
end

function assertNumericalRank(matrix, expectedRank, message)
% Сингулярные числа отделяют физические формы от численных нулевых форм.
singularValues = svd(matrix);
tolerance = max(singularValues) * 1e-10;
actualRank = sum(singularValues > tolerance);
if actualRank ~= expectedRank
    error('MKEF:VerificationFailed', '%s Rank=%d, expected=%d.', ...
        message, actualRank, expectedRank);
end
end

function assertThrows(expectedIdentifier, operation)
% Требуем ожидаемую ошибку проверки данных, а не произвольный сбой.
try
    operation();
catch exception
    if strcmp(exception.identifier, expectedIdentifier)
        return;
    end
    error('MKEF:VerificationFailed', ...
        'Expected error %s, received %s.', ...
        expectedIdentifier, exception.identifier);
end
error('MKEF:VerificationFailed', 'Expected error %s.', expectedIdentifier);
end

function assertClose(actual, expected, relativeTolerance, absoluteTolerance, message)
% Сочетание допусков подходит и для нулевых элементов, и для больших жёсткостей.
difference = max(abs(actual(:) - expected(:)));
scale = max(abs(expected(:)));
limit = absoluteTolerance + relativeTolerance * scale;
if difference > limit
    error('MKEF:VerificationFailed', '%s Error=%g, tolerance=%g.', ...
        message, difference, limit);
end
end
