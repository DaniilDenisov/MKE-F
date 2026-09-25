% Базовый класс структурного КЭ.
%
% Octave не поддерживает используемый MATLAB-синтаксис прототипов
% абстрактных методов вне @-каталогов. Поэтому базовые реализации явно
% сообщают об ошибке, а классы элементов переопределяют эти методы.
classdef FiniteElementStructural < handle
    properties (Access = public)
        % Тип элемента.
        elType;
    end
    properties (Access = protected)
        % Координаты узлов.
        elNodesCoords;
        % Номера узлов.
        elNodesNums;
        % Данные элемента (плотность, площадь).
        elData;
    end
    methods (Access = public)
        % Конструктор без аргументов.
        function obj = FiniteElementStructural()
        end
        % Функция получения всех узловых координат (для печати сетки).
        function nCoords = GetNodalCoords(this)
            nCoords = this.elNodesCoords;
        end
        % Функция получения номеров узлов элемента (для построения матрицы
        % соответствия).
        function nNums = GetNodesNums(this)
            nNums = this.elNodesNums;
        end
    end
    methods (Access = public)
        % Ассемблер эл-та в глоб. матрицы M и K по матрице соответствия IM.
        function [GK, GM] = Assembler(~, GK, GM, ~)
            error('FiniteElementStructural:NotImplemented', ...
                'Assembler must be implemented by an element subclass.');
        end
    end
    methods (Access = protected)
        % Вычисление матриц массы и жесткости.
        function M = MassElementMatrix(~)
            M = [];
            error('FiniteElementStructural:NotImplemented', ...
                'MassElementMatrix must be implemented by an element subclass.');
        end
        function K = StiffnessElementMatrix(~)
            K = [];
            error('FiniteElementStructural:NotImplemented', ...
                'StiffnessElementMatrix must be implemented by an element subclass.');
        end
    end
end
