% Класс КЭ сетки.
% Copyright 2017 Daniil S. Denisov
classdef FEMesh < handle
    properties (Access = public)
        % Число элементов.
        numberOfElems=0;
        % Число узлов.
        numberOfNodes=0;
        % Число ГУ с приложенной силой.
        numberOfForceBCs=0;
        % Число гу заделки.
        numberOfFixBCs=0;
        % Число степеней свободы на один узел.
        dofPerNode;
        % Поузловая матрица соответствия.
        iMnod;
        % Массив элементов (пользовательский класс).
        allMeshElems;
        % Массив узлов.
        allNodes;
        % Массив ГУ заделки.
        allFixBCs;
        % Массив ГУ приложенных сил.
        allForceBCs;
    end
    methods
        % Чтение файла сетки (конструктор).
        function obj = FEMesh(filename)
            % Создание дескриптора файла. -1 если файла нет.
            fid = fopen(filename,'r');
            % Макс. число блоков = 10.
            for bn=1:10
                % Чтение строки.
                line = fgetl(fid);
                % Первая строка должна быть маркером. Решение по типу маркера.
                switch line
                    case 'nodes'
                        obj.readNodes(fid);
                    case 'elems_112'
                        obj.readElems112(fid);
                    case 'elems_113'
                        obj.readElems113(fid);
                    case 'bcforce_stat'
                        obj.readBcForcesStat(fid);
                    case 'bcforce_harm'
                        obj.readBcForcesHarm(fid);
                    case 'bcfix'
                        obj.readBcFix(fid);
                    case -1
                        break
                    otherwise
                        warning('Unexpected line marker in task file.');
                end
            end
            % Закрытие файла.
            fclose(fid);
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
            disp(this.iMnod);
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
        % Печать элементов.
        function DispElems(this)
            for i=1:this.numberOfElems
               this.allMeshElems(i).Disp();
            end
        end
        
        % -------------------- Блоки данных (чтение) --------------------
        % Чтение блока данных об узлах.
        function readNodes(this,fid)
            % В начале блока должен быть размер.
            line = fgetl(fid);
            nNum = sscanf(line,'%d');
            this.numberOfNodes = nNum;
            % Чтение строк и сохранение поле узлов.
            for i=1:nNum
                this.allNodes(i,:) = sscanf(fgetl(fid),'%f,%f,%f');
            end
        end
        
        % Чтение блока данных об элементах 112.
        % Внимание! Массив узлов должен быть заполнен.
        function readElems112(this, fid)
            % В начале блока должен быть размер.
            line = fgetl(fid);
            elNum = sscanf(line,'%d');
            this.numberOfElems = elNum;
            % Преаллокация набора пустых КЭ типа 112 (Truss2DElement).
            elemArrTmp(elNum,1) = Truss2DElement();
            this.allMeshElems = elemArrTmp;
            % Для элемента 112 должно быть 2 СС в узле.
            this.dofPerNode = 2;
            % Далее - данные. В соответствии с размером.
            for i=1:elNum
                line = fgetl(fid);
                splitElemLine = strsplit(line,',');
                % Чтение типа элемента пропущено, т.к. реализованы маркеры.
                % Чтение со второй позиции. Номера узлов.
                elNodes = str2double(splitElemLine(2:3));
                elData = str2double(splitElemLine(4:6));
                elNode1Coords = this.allNodes(elNodes(1),:);
                elNode2Coords = this.allNodes(elNodes(2),:);
                this.allMeshElems(i).SetupElement(...
                            [elNode1Coords; elNode2Coords],...
                            elNodes, ...
                            elData);
            end
            % Можем создать матрицу соответствия, т.к. есть dofPerNode.
            for i=1:this.numberOfNodes
                this.iMnod(i,:) = [(i*this.dofPerNode)-1 i*this.dofPerNode];
            end
        end
        
        % Чтение блока данных об элементах 113 (балка 2 узла по 3 СС).
        % Внимание! Массив узлов должен быть заполнен.
        function readElems113(this, fid)
            % В начале блока должен быть размер.
            line = fgetl(fid);
            elNum = sscanf(line,'%d');
            this.numberOfElems = elNum;
            % Преаллокация набора пустых КЭ типа 112 (Truss2DElement).
            elemArrTmp(elNum,1) = Beam2DElement();
            this.allMeshElems = elemArrTmp;
            % Для элемента 113 должно быть 3 СС в узле.
            this.dofPerNode = 3;
            % Далее - данные. В соответствии с размером.
            for i=1:elNum
                line = fgetl(fid);
                splitElemLine = strsplit(line,',');
                % Чтение типа элемента пропущено, т.к. реализованы маркеры.
                % Чтение со второй позиции. Номера узлов.
                elNodes = str2double(splitElemLine(2:3));
                % Данные элемента балки (площадь, мод. Юнга, плотность,
                % момент инерции)
                elData = str2double(splitElemLine(4:7));
                elNode1Coords = this.allNodes(elNodes(1),:);
                elNode2Coords = this.allNodes(elNodes(2),:);
                this.allMeshElems(i).SetupElement(...
                            [elNode1Coords; elNode2Coords],...
                            elNodes, ...
                            elData);
            end
            % Можем создать матрицу соответствия, т.к. есть dofPerNode.
            for i=1:this.numberOfNodes
                this.iMnod(i,:) = [(i*this.dofPerNode)-2 ...
                                   (i*this.dofPerNode)-1 ...
                                   (i*this.dofPerNode)];
            end
        end
        
        % Чтение блока данных о граничных условиях заделки.
        function readBcFix(this, fid)
            % В начале блока должен быть размер.
            line = fgetl(fid);
            bcFixNum = sscanf(line,'%d');
            this.numberOfFixBCs = this.numberOfFixBCs + bcFixNum;
            % Чтение строк и сохранение поле узлов.
            for i=1:bcFixNum
                this.allFixBCs(i,:) = sscanf(fgetl(fid),...
                    '%f,%f,%f,%f,%f');
            end
        end
        % Чтение блока данных о ГУ прилож. сил (стат./удар).
        function readBcForcesStat(this, fid)
            % В начале блока должен быть размер.
            line = fgetl(fid);
            bcForceNum = sscanf(line,'%d');
            this.numberOfForceBCs = this.numberOfForceBCs + bcForceNum;
            % Чтение строк и сохранение поле узлов.
            for i=1:bcForceNum
                this.allForceBCs(i,:) = sscanf(fgetl(fid),...
                    '%f,%f,%f,%f,%f');
            end
        end
        % Чтение блока данных о ГУ приложенных сил (гармоническая).
        function readBcForcesHarm(this, fid)
            % В начале блока должен быть размер.
            line = fgetl(fid);
            bcForceNum = sscanf(line,'%d');
            this.numberOfForceBCs = this.numberOfForceBCs + bcForceNum;
            % Чтение строк и сохранение поле узлов.
            for i=1:bcForceNum
                this.allForceBCs(i,:) = sscanf(fgetl(fid),...
                    '%f,%f,%f,%f,%f,%f');
            end
        end
    end
end

