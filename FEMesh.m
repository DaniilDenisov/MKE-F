% Класс КЭ сетки.
% Copyright 2017 Daniil S. Denisov
classdef FEMesh < handle
    properties (Access = public)
        numberOfElems;
        numberOfNodes;
        numberOfBCs;
        dofPerNode;
        % Поэлементная м-ца соответствия. Узловые СС - к глоб номерам СС.
        iM;
        % Поузловая матрица соответствия.
        iMnod;
        allMeshElems;
        allNodes;
        allBCs;
    end
    methods
        % Чтение файла сетки (конструктор).
        function obj = FEMesh(filename)
            % Чтение строки с размерами массивов.
            totalLine = dlmread(filename,',',[0 0 0 2]);
            % Разбор строки с размерами. Вывод количеств "всего".
            obj.numberOfElems = totalLine(1);
            obj.numberOfNodes = totalLine(2);
            obj.numberOfBCs = totalLine(3);
            fprintf('Els:%d Nodes:%d BCs:%d\n', obj.numberOfElems,...
                obj.numberOfNodes, obj.numberOfBCs);
            % Чтение строк данных в соответствии с размерами массивов.
            % Тип эл-та, номера узлов в элементах и данные элементов.
            allElemsTmp = dlmread(filename,',',...
                [1 0 obj.numberOfElems 6]);
            % Координаты узлов.
            obj.allNodes = dlmread(filename,',',...
                [obj.numberOfElems+1 0 ...
                obj.numberOfNodes+obj.numberOfElems 2]);
            % Граничные условия.
            obj.allBCs = dlmread(filename,',',...
                [obj.numberOfNodes+obj.numberOfElems+1 0 ...
                obj.numberOfNodes+obj.numberOfElems+obj.numberOfBCs 5]);
            % Преаллокация полного набора пустых КЭ.
            elemArrTmp(obj.numberOfElems,1) = Truss2DElement();
            obj.allMeshElems = elemArrTmp;
            % Наполнение пустых КЭ данными.
            for i=1:obj.numberOfElems
                % Выбрать строку с данными i-го элемента.
                currElem = allElemsTmp(i,:);
                % Взять тип элемента.
                currType = currElem(1);
                % Взять номера узлов. Два номера узлов.
                currNode1 = currElem(2);
                currNode2 = currElem(3);
                % Поля данных: площадь, модуль Юнга, плотность,
                % момент инерции поперечного сечения.
                currArea = currElem(4);
                currEmod = currElem(5);
                currRho = currElem(6);
                currI = currElem(7);   % Для элемента фермы =0.
                % Координаты узлов.
                currNode1_XYZ = [obj.allNodes(currNode1,1) ...
                    obj.allNodes(currNode1,2) ...
                    obj.allNodes(currNode1,3)];
                currNode2_XYZ = [obj.allNodes(currNode2,1) ...
                    obj.allNodes(currNode2,2) ...
                    obj.allNodes(currNode2,3)];
                % Наполнение объекта элемента в объекте FEМеsh и
                % заполнение строки таблицы соответствия iM.
                switch currType
                    % Для элемента типа 112.
                    case 112
                        obj.allMeshElems(i).SetupElement(...
                            [currNode1_XYZ; currNode2_XYZ],...
                            [currNode1,currNode2], ...
                            [currArea, currEmod, currRho]);
                        % Для матрицы соответствия эл-та тип 112.
                        obj.dofPerNode = 2;
                        % Взять номера узлов.
                        nodeNums = obj.allMeshElems(i).GetNodesNums;
                        % Вычислить глобальные индексы степеней свободы.
                        GD1 = (nodeNums(1)-1)*(obj.dofPerNode)+1;
                        GD2 = GD1+1;
                        GD3 = (nodeNums(2)-1)*(obj.dofPerNode)+1;
                        GD4 = GD3+1;
                        obj.iM(i,1) = GD1;
                        obj.iM(i,2) = GD2;
                        obj.iM(i,3) = GD3;
                        obj.iM(i,4) = GD4;
                        % Создание (перезапись) записи в поузловой матрице
                        % соответствия iMnod. Не эффективно.
                        % Присваиваются повторно. Надо доработать.
                        obj.iMnod(nodeNums(1),1) = GD1;
                        obj.iMnod(nodeNums(1),2) = GD2;
                        obj.iMnod(nodeNums(2),1) = GD3;
                        obj.iMnod(nodeNums(2),2) = GD4;
                    case 113
                        obj.allMeshElems(i).SetupElement(...
                            [currNode1_XYZ; currNode2_XYZ],...
                            [currNode1,currNode2], ...
                            [currArea, currEmod, currRho, currI]);
                        % Для матрицы соответствия эл-та тип 112.
                        obj.dofPerNode = 3;
                        % Взять номера узлов.
                        nodeNums = obj.allMeshElems(i).GetNodesNums;
                        % Вычислить глобальные индексы степеней свободы.
                        GD1 = (nodeNums(1)-1)*(obj.dofPerNode)+1;
                        GD2 = GD1+1;
                        GD3 = GD2+1;
                        GD4 = (nodeNums(2)-1)*(obj.dofPerNode)+1;
                        GD5 = GD4+1;
                        GD6 = GD5+1;
                        % Заполнение элементной матрицы соответствия.
                        obj.iM(i,1) = GD1;
                        obj.iM(i,2) = GD2;
                        obj.iM(i,3) = GD3;
                        obj.iM(i,4) = GD4;
                        obj.iM(i,5) = GD5;
                        obj.iM(i,6) = GD6;
                        % Создание (перезапись) записи в поузловой матрице
                        % соответствия iMnod. Не эффективно.
                        % Если узел уже был в другом эл-те, то
                        % тут присваиваются повторно. Надо доработать.
                        obj.iMnod(nodeNums(1),1) = GD1;
                        obj.iMnod(nodeNums(1),2) = GD2;
                        obj.iMnod(nodeNums(1),3) = GD3;
                        obj.iMnod(nodeNums(2),1) = GD4;
                        obj.iMnod(nodeNums(2),2) = GD5;
                        obj.iMnod(nodeNums(2),3) = GD6;
                    otherwise
                        disp('FEMesh: Unknown element type.');
                end
            end
        end
        % Изображение 2D сетки на графике.
        function Plot2DMesh(this)
            labels = cellstr(num2str([1:this.numberOfNodes]'));
            plot(this.allNodes(:,1),this.allNodes(:,2), 'bo',...
                'MarkerSize',8,'MarkerFaceColor', 'w',...
                'LineWidth', 2.0);
            % Поля на графике:
            nodesSpanX = max(this.allNodes(:,1))-min(this.allNodes(:,1));
            nodesSpanY = max(this.allNodes(:,2))-min(this.allNodes(:,2));
            % Если не все элементы в ряд.
            if (nodesSpanX~=0 && nodesSpanY~=0)
                xlim([(-1)*((nodesSpanX)/10) ...
                    max(this.allNodes(:,1))+(nodesSpanX)/10]);
                ylim([(-1)*((nodesSpanY)/10) ...
                    max(this.allNodes(:,2))+(nodesSpanY)/10]);
            end
            % Печать марок.
            text(this.allNodes(:,1),this.allNodes(:,2), labels,...
                'LineWidth',2,'VerticalAlignment', 'bottom',...
                'HorizontalAlignment','right');
            % Печать элементов.
            hold on;
            grid on;
            for el = 1:this.numberOfElems
                % Определение координат узлов текущего элемента.
                currElem = this.allMeshElems(el);
                currNodalCoords = currElem.GetNodalCoords();
                % Печать линии.
                plot([currNodalCoords(1,1) currNodalCoords(2,1)],...
                    [currNodalCoords(1,2) currNodalCoords(2,2)]);
                % Определение координат текста обозначения элемента.
                dx = currNodalCoords(1,1) - currNodalCoords(2,1);
                dy = currNodalCoords(1,2) - currNodalCoords(2,2);
                textCoord1 = currNodalCoords(1,1)-(dx/2);
                textCoord2 = currNodalCoords(1,2)-(dy/2);
                % Внесение на график обозначения элемента
                text(textCoord1,textCoord2,num2str(el));
            end
            hold off;
        end
        % Печать матрицы соответствия.
        function DispIM(this)
            format shortG;
            disp('DOF distribution matrix:');
            disp(this.iM);
            format compact;
        end
        % Печать узлов из элементов сетки.
        function DispNN(this)
            format shortG;
            disp('Nodes:');
            for i=1:this.numberOfElems
                nn = this.allMeshElems(i).GetNodesNums();
                disp(nn);
            end
            format compact;
        end
    end
end

