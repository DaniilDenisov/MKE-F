% Абстрактный класс структурного КЭ.
classdef (Abstract) FiniteElementStructural < handle
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
    methods (Access = public, Abstract = true)
        % Ассемблер эл-та в глоб. матрицы M и K по матрице соответствия IM.
        [GM, GK] = Assembler(GM, GK, IM, elNum)
    end
    methods (Access = protected, Abstract = true)
        % Вычисление матриц массы и жесткости.
        M = MassElementMatrix(this)
        K = StiffnessElementMatrix(this)
    end
end