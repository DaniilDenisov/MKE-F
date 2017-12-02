% Класс балочного элемента плоской рамы.
% Copyright 2017 Daniil S. Denisov
classdef Beam2DElement < FiniteElementStructural
    methods (Access = public)
        % Реализация абстрактного метода ансамблирования.
        function [GK, GM] = Assembler(this, GK, GM, IM, elNum)
%             % Вычисление элементных МЖ и ММ.
%             K = StiffnessElementMatrix(this);
%             M = MassElementMatrix(this);
%             % Опред. глоб. индексов для элементных степеней свободы.
%             glDOF1 = IM(elNum,1);
%             glDOF2 = IM(elNum,2);
%             glDOF3 = IM(elNum,3);
%             glDOF4 = IM(elNum,4);
%             % Ансамблирование элементов ke в глоб. матрицу GK.
%             GK(glDOF1,glDOF1) = GK(glDOF1,glDOF1) + K(1,1);
%             GK(glDOF1,glDOF2) = GK(glDOF1,glDOF2) + K(1,2);
%             GK(glDOF1,glDOF3) = GK(glDOF1,glDOF3) + K(1,3);
%             GK(glDOF1,glDOF4) = GK(glDOF1,glDOF4) + K(1,4);
%             GK(glDOF2,glDOF1) = GK(glDOF2,glDOF1) + K(2,1);
%             GK(glDOF2,glDOF2) = GK(glDOF2,glDOF2) + K(2,2);
%             GK(glDOF2,glDOF3) = GK(glDOF2,glDOF3) + K(2,3);
%             GK(glDOF2,glDOF4) = GK(glDOF2,glDOF4) + K(2,4);
%             GK(glDOF3,glDOF1) = GK(glDOF3,glDOF1) + K(3,1);
%             GK(glDOF3,glDOF2) = GK(glDOF3,glDOF2) + K(3,2);
%             GK(glDOF3,glDOF3) = GK(glDOF3,glDOF3) + K(3,3);
%             GK(glDOF3,glDOF4) = GK(glDOF3,glDOF4) + K(3,4);
%             GK(glDOF4,glDOF1) = GK(glDOF4,glDOF1) + K(4,1);
%             GK(glDOF4,glDOF2) = GK(glDOF4,glDOF2) + K(4,2);
%             GK(glDOF4,glDOF3) = GK(glDOF4,glDOF3) + K(4,3);
%             GK(glDOF4,glDOF4) = GK(glDOF4,glDOF4) + K(4,4);
%             % Ансамблирование матрицы жесткости в глобальную GM.
%             GM(glDOF1,glDOF1) = GM(glDOF1,glDOF1) + M(1,1);
%             GM(glDOF1,glDOF2) = GM(glDOF1,glDOF2) + M(1,2);
%             GM(glDOF1,glDOF3) = GM(glDOF1,glDOF3) + M(1,3);
%             GM(glDOF1,glDOF4) = GM(glDOF1,glDOF4) + M(1,4);
%             GM(glDOF2,glDOF1) = GM(glDOF2,glDOF1) + M(2,1);
%             GM(glDOF2,glDOF2) = GM(glDOF2,glDOF2) + M(2,2);
%             GM(glDOF2,glDOF3) = GM(glDOF2,glDOF3) + M(2,3);
%             GM(glDOF2,glDOF4) = GM(glDOF2,glDOF4) + M(2,4);
%             GM(glDOF3,glDOF1) = GM(glDOF3,glDOF1) + M(3,1);
%             GM(glDOF3,glDOF2) = GM(glDOF3,glDOF2) + M(3,2);
%             GM(glDOF3,glDOF3) = GM(glDOF3,glDOF3) + M(3,3);
%             GM(glDOF3,glDOF4) = GM(glDOF3,glDOF4) + M(3,4);
%             GM(glDOF4,glDOF1) = GM(glDOF4,glDOF1) + M(4,1);
%             GM(glDOF4,glDOF2) = GM(glDOF4,glDOF2) + M(4,2);
%             GM(glDOF4,glDOF3) = GM(glDOF4,glDOF3) + M(4,3);
%             GM(glDOF4,glDOF4) = GM(glDOF4,glDOF4) + M(4,4);
        end
        % Конструктор класса без аргументов.
        function obj = Beam2DElement()
        end
        % Функция установки полей.
        function SetupElement(this,nodCoordsIn,...
                nodesNumsIn ,dataIn)
            % Установка типа элемента для конструктора.
            this.elType = 112;
            this.elNodesCoords = nodCoordsIn;
            this.elData = dataIn;
            this.elNodesNums = nodesNumsIn;
        end
        % Функция печати.
        function Disp(this)
            format shortG;
            fprintf('type:%d\n',this.elType);
            disp(this.elNodesCoords);
            disp(this.elNodesNums);
            disp(this.elData);
            format compact;
        end
    end
    methods (Access = protected)
        % Функция определения матрицы масс элемента.
        function M = MassElementMatrix(this)
%             % Получение из "поля данных" характеристик элемента.
%             currArea = this.elData(1);
%             currRho = this.elData(3);
%             % Вызов функции определения матрицы косинусов и длины.
%             [T, length] = TransformMatrix(this);
%             % Матрица масс элемента без преобразования координат.
%             meInit = (currRho*currArea*length/6)*...
%                 [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
%             % Преобразование матрицы масс.
%             M = T'*meInit*T;
        end
        % Функция определения матрицы жесткости элемента.
        function K = StiffnessElementMatrix(this)
%             % Получение из "поля данных" характеристик элемента.
%             currArea = this.elData(1);
%             currEmod = this.elData(2);
%             % Вызов функции определения матрицы косинусов и длины.
%             [T, length] = TransformMatrix(this);
%             % Матрица жесткости элемента без преобразования координат.
%             KInit = zeros(4,4);
%             kCoeff = currArea*currEmod/length;
%             KInit(1,1) = kCoeff;
%             KInit(1,3) = -kCoeff;
%             KInit(3,1) = -kCoeff;
%             KInit(3,3) = kCoeff;
%             % Преобразование элементной матрицы жесткости.
%             K = T'*KInit*T;
        end
        % Функция определения матрицы косинусов и длины элемента.
        function [T, length] = TransformMatrix(this)
            node1 = this.elNodesCoords(1,:);
            node2 = this.elNodesCoords(2,:);
            dx = node2(1,1)-node1(1,1);
            dy = node2(1,2)-node1(1,2);
            length = sqrt(dx^2 + dy^2);
            c = dx/length;
            s = dy/length;
            T = [c s 0 0;-s c 0 0;0 0 c s;0 0 -s c];
        end
    end
end