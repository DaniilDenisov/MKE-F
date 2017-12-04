% Класс задачи.
% Copyright 2017 Daniil S. Denisov
classdef StructFEProblem < handle
    properties (Access = public)
        % Имя кейс-файла.
        filename
        % Матрица жесткости.
        K
        % Матрица масс.
        M
        % Вектор правой части.
        F
        % Сетка с ГУ.
        mesh
    end
    methods
        % Конструктор с аргументом.
        function obj = StructFEProblem(filename)
            obj.filename = filename;
            % Чтение кейса, создание сетки.
            obj.mesh = FEMesh(obj.filename);
            % Вывод сетки на график.
            obj.mesh.Plot2DMesh();
            % Печать номеров узлов в элементах.
            obj.mesh.DispNN();
            % Печать матрицы соответствия.
            obj.mesh.DispIM();
            % Создание глобальной матрицы жесткости, вектора F и матрицы
            % масс в зависимости от типа элемента.
            elType = obj.mesh.allMeshElems(1).elType;
            switch elType
                % Элемент Тип 112:
                case 112
                    % Степеней свободы в узле.
                    dofPerNode = 2;
                    % Определение общего числа степеней свободы в системе
                    % и иниц. глобальной матрицы жесткости.
                    systemDOF = obj.mesh.numberOfNodes*dofPerNode;
                    GlobK = zeros(systemDOF,systemDOF);
                    GlobM = zeros(systemDOF,systemDOF);
                    % Обход всех элементов в сетке.
                    for i=1:obj.mesh.numberOfElems
                        % Ансамблирование в глобальные матрицы.
                        [GlobK,GlobM] = ...
                            obj.mesh.allMeshElems(i).Assembler(GlobK,...
                            GlobM, obj.mesh.iM, i);
                    end
                    obj.K = GlobK;
                    obj.M = GlobM;
                    % Вектор правой части (сил). Преаллокация без ГУ.
                    obj.F = zeros(obj.mesh.numberOfNodes*dofPerNode,1);
                % Элемент Тип 113.
                case 113
                    % Степеней свободы в узле.
                    dofPerNode = 3;
                    % Определение общего числа степеней свободы в системе
                    % и иниц. глобальной матрицы жесткости.
                    systemDOF = obj.mesh.numberOfNodes*dofPerNode;
                    GlobK = zeros(systemDOF,systemDOF);
                    GlobM = zeros(systemDOF,systemDOF);
                    % Обход всех элементов в сетке.
                    for i=1:obj.mesh.numberOfElems
                        % Ансамблирование в глобальные матрицы.
                        [GlobK,GlobM] = ...
                            obj.mesh.allMeshElems(i).Assembler(GlobK,...
                            GlobM, obj.mesh.iM, i);
                    end
                    obj.K = GlobK;
                    obj.M = GlobM;
                    % Вектор правой части (сил). Преаллокация без ГУ.
                    obj.F = zeros(obj.mesh.numberOfNodes*dofPerNode,1);
                otherwise
                    disp('SructFEProblem(filename): bad element type!');
            end
            
        end
        % Метод наложения ГУ. Для статики без матрицы масс.
        function ApplyBC(this)
            BC = this.mesh.allBCs;
            BCtotal = this.mesh.numberOfBCs;
            for i=1:BCtotal
                % Взять ГУ
                currBC = BC(i,:);
                % Взять тип ГУ.
                typeBC = currBC(1);
                % Взять номер узла.
                nnumBC = currBC(2);
                % Вычислить номера степеней свободы в глоб. МЖ, MM, ВПЧ.
                GLDOFs = this.mesh.iMnod(nnumBC,:);
                % Применить ГУ к МЖ.
                % Заделка всех СС. Обнулить строки, столбцы и записать 1 на
                % диагональ.
                if (typeBC==1)
                    for n=1:size(GLDOFs,1)
                        this.K(GLDOFs(n),:) = 0;
                        this.K(:,GLDOFs(n)) = 0;
                        this.K(GLDOFs(n),GLDOFs(n)) = 1;
                    end
                end
                % Заделка вертик. СС, обнулить один (по второй СС) столбец.
                if (typeBC==2)
                    this.K(GLDOFs(2),:) = 0;
                    this.K(:,GLDOFs(2)) = 0;
                    this.K(GLDOFs(2),GLDOFs(2)) = 1;
                end
                % Заделка гор. СС, обнулить один (по первой СС) столбец.
                if (typeBC==3)
                    this.K(GLDOFs(1),:) = 0;
                    this.K(:,GLDOFs(1)) = 0;
                    this.K(GLDOFs(1),GLDOFs(1)) = 1;
                end
                % Если сила, то установить значение в вектор правой части.
                if (typeBC==10)
                    % Взять величину силы по компонентам из BC.
                    forceValue = zeros(3,1);
                    forceValue(1) = currBC(1,3);
                    forceValue(2) = currBC(1,4);
                    forceValue(3) = currBC(1,5);
                    % Установить в вектор пр. части.
                    for n=1:size(GLDOFs,1)
                        this.F(GLDOFs(1)) = this.F(GLDOFs(1))+forceValue(1);
                        this.F(GLDOFs(2)) = this.F(GLDOFs(2))+forceValue(2);
                    end
                end
            end
        end
        % Метод наложения ГУ для расчетов собственных колебаний. Без ВПЧ.
        function ApplyBCEig(this)
            BC = this.mesh.allBCs;
            BCtotal = this.mesh.numberOfBCs;
            for i=1:BCtotal
                % Взять ГУ
                currBC = BC(i,:);
                % Взять тип ГУ.
                typeBC = currBC(1);
                % Взять номер узла.
                nnumBC = currBC(2);
                % Вычислить номер степени свободы в глоб. МЖ, MM, ВПЧ.
                GLDOFs = this.mesh.iMnod(nnumBC,:);
                glDOF1=GLDOFs(1);
                glDOF2=GLDOFs(2);
                % Применить ГУ к МЖ.
                % Заделка обоих СС. Обнулить строки, столбцы и записать 1 на диагональ.
                % Т.к. в узле две СС, то обнуляется 2 строки и 2 столбца.
                if (typeBC==1)
                    this.M(glDOF1,:) = 0;
                    this.M(:,glDOF1) = 0;
                    this.M(glDOF2,:) = 0;
                    this.M(:,glDOF2) = 0;
                    this.M(glDOF1,glDOF1) = 1;
                    this.M(glDOF2,glDOF2) = 1;
                    this.K(glDOF1,:) = 0;
                    this.K(:,glDOF1) = 0;
                    this.K(glDOF2,:) = 0;
                    this.K(:,glDOF2) = 0;
                end
                % Заделка вертик. СС, обнулить один (по второй СС) столбец.
                if (typeBC==2)
                    this.M(glDOF2,:) = 0;
                    this.M(:,glDOF2) = 0;
                    this.M(glDOF2,glDOF2) = 1;
                    this.K(glDOF2,:) = 0;
                    this.K(:,glDOF2) = 0;
                end
                % Заделка гор. СС, обнулить один (по первой СС) столбец.
                if (typeBC==3)
                    this.M(glDOF1,:) = 0;
                    this.M(:,glDOF1) = 0;
                    this.M(glDOF1,glDOF1) = 1;
                    this.K(glDOF1,:) = 0;
                    this.K(:,glDOF1) = 0;
                end
            end
        end
        % Метод наложения ГУ для динамического расчета (Transient).
        function ApplyBCTrt(this,tStep,tsNum)
            % this  - объект для работы с полями (неявный аргумент),
            % tStep - шаг по времени,
            % tsNum - количество временных шагов.
            % Копирование в локальные переменные числа узлов и СС на узел.
            nn = this.mesh.numberOfNodes;
            dpn = this.mesh.dofPerNode;
            % Инициализация вектора правой части для динамического расчета.
            % Столбец - шаг по времени.
            this.F = zeros(nn*dpn,tsNum);
            % Наложение ограничений.
            allBCs = this.mesh.allBCs;
            for bc=1:this.mesh.numberOfBCs
                % Взять ГУ
                currBC = allBCs(bc,:);
                % Взять тип ГУ.
                typeBC = currBC(1);
                % Взять номер узла.
                nnumBC = currBC(2);
                % Вычислить номер степени свободы в глоб. МЖ, MM, ВПЧ.
                GLDOFs = this.mesh.iMnod(nnumBC,:);
                glDOF1=GLDOFs(1);
                glDOF2=GLDOFs(2);
                % Наложение ГУ различных типов.
                % ГУ перемещения на узел.
                % Заделка обоих СС. Обнулить строки, столбцы и записать 1 на диагональ.
                % Т.к. в узле две СС, то обнуляется 2 строки и 2 столбца.
                if (typeBC==1)
                    this.K(glDOF1,:) = 0;
                    this.K(:,glDOF1) = 0;
                    this.K(glDOF2,:) = 0;
                    this.K(:,glDOF2) = 0;
                    this.K(glDOF1,glDOF1) = 1;
                    this.K(glDOF2,glDOF2) = 1;
                end
                % Заделка вертик. СС, обнулить один (по второй СС) столбец.
                if (typeBC==2)
                    this.K(glDOF2,:) = 0;
                    this.K(:,glDOF2) = 0;
                    this.K(glDOF2,glDOF2) = 1;
                end
                % Заделка гор. СС, обнулить один (по первой СС) столбец.
                if (typeBC==3)
                    this.K(glDOF1,:) = 0;
                    this.K(:,glDOF1) = 0;
                    this.K(glDOF1,glDOF1) = 1;
                end
                % Если сила, то установить значение в вектор правой части
                % на каждом шаге по времени (каждом столбце).
                if (typeBC==10)
                    % Взять величину силы по компонентам из BC.
                    forceValue = zeros(3,1);
                    forceValue(1) = currBC(1,3);
                    forceValue(2) = currBC(1,4);
                    forceValue(3) = currBC(1,5);
                    % Установить в вектор пр. части.
                    this.F(glDOF1,:) = this.F(glDOF1,:)+forceValue(1);
                    this.F(glDOF2,:) = this.F(glDOF2,:)+forceValue(2);
                end
                % Если прикладывается сила гармоническая, создать tsNum
                % столбцов с сохранением неизменных сил.
                if typeBC==11
                    forceValue = zeros(3,1);
                    % Считать Fx,Fy,Fz в forceValue.
                    forceValue(1) = currBC(1,3);
                    forceValue(2) = currBC(1,4);
                    forceValue(3) = currBC(1,5);
                    % Считать частоту в forceValue.
                    freq = currBC(1,6);
                    % Формирование вектора правой части для каждого шага
                    % по времени.
                    for i=1:tsNum
                        this.F(glDOF1,i) = this.F(glDOF1,i)+...
                            forceValue(1)*sin((2*pi*freq)*(i*tStep));
                        this.F(glDOF2,i) = this.F(glDOF2,i)+...
                            forceValue(2)*sin((2*pi*freq)*(i*tStep));
                    end
                end
            end
        end
        % Метод запуска расчета статического нагружения.
        function RunStatic(this)
            % --> Расчет статики.
            % Сохранение МЖ и ММ без ГУ. Например, для вычисления реакций.
            KnoBC = this.K;
            MnoBC = this.M;
            % Наложение ГУ.
            this.ApplyBC();
            % Решение системы, определение перемещений. Вывод.
            dspl = this.K\this.F;
            disp('DOFs (displ. components):')
            disp(dspl);
            % Определение реакций. Вывод.
            reacts = KnoBC*dspl-this.F;
            disp('Reactions (components):')
            disp(reacts);
            % Восстановление МЖ и ММ для следующих расчетов.
            this.K = KnoBC;
            this.M = MnoBC;
        end
        % Метод запуска расчета собственных колебаний.
        function RunModal(this)
            % --> Расчет собственных колебаний.
            % Сохранение МЖ и ММ без ГУ.
            KnoBC = this.K;
            MnoBC = this.M;
            % Наложение ГУ.
            this.ApplyBCEig();
            % Решение модальной задачи.
            freqSol = (1/(2*pi))*sqrt(eig(this.K, this.M));
            disp('Natural frequencies (Hz):')
            disp(freqSol);
            % Восстановление МЖ и ММ для следующих расчетов.
            this.K = KnoBC;
            this.M = MnoBC;
        end
        % Метод запуска анализа динамики со внешними силами.
        function RunTransient(this,tStep,tDur,node,dofToPlot)
            % this  - объект для работы с полями (неявный аргумент),
            % tStep - шаг по времени,
            % tDur  - время моделирования общее,
            % node  - узел в котором строится график по времени.
            % doftoplot - номер СС в узле для вывода.
            % Целое число временных шагов.
            tsNum = fix(tDur/tStep);
            % Массивы для узловых скоростей, ускорений и перемещений.
            nn = this.mesh.numberOfNodes;
            dpn = this.mesh.dofPerNode;
            acsNodal = zeros(nn*dpn,1);
            speNodal = zeros(nn*dpn,1);
            dspNodal = zeros(nn*dpn,1);
            % Установка вектора правой части в ноль для всех шагов.
            this.F = zeros(nn*dpn,tsNum);
            % ГУ.
            this.ApplyBCTrt(tStep,tsNum);
            % Постоянные для прямого интегрирования Ньюмарка.
            delta = 0.5;
            alfa = 0.25;
            a0 = 1/(alfa*(tStep^2));
            a1 = delta/(alfa*(tStep));
            a2 = 1/(alfa*(tStep));
            a3 = (1/(2*alfa))-1;
            a4 = (delta/alfa)-1;
            a5 = (tStep/2)*((delta/alfa)-2);
            a6 = tStep*(1-delta);
            a7 = delta*tStep;
            % Наложение ГУ на МЖ и ММ.
            this.ApplyBC();
            % Формирование эффективной матрицы жесткости.
            Khat = this.K+a0*this.M;
            % Для каждого шага по времени считаем перемещения.
            for i=1:tsNum
                % Эффективная нагрузка.
                Rhat = this.F(:,i)+this.M*(a0*dspNodal(:,i)+a2*speNodal(:,i)+...
                    a3*acsNodal(:,i));
                % Определяем перемещения для следующего шага по времени.
                dspNodal(:,i+1) = Khat\Rhat;
                % Определить ускорения.
                acsNodal(:,i+1) = a0*(dspNodal(:,i+1)-dspNodal(:,i))-...
                    a2*speNodal(:,i)-a3*acsNodal(:,i);
                % Определить скорости.
                speNodal(:,i+1) = speNodal(:,i)+a6*acsNodal(:,i)+...
                    a7*acsNodal(:,i+1);
            end
            % Переменная времени.
            timeSpan = 0:tStep:tsNum*tStep;
            % Печать перемещений.
            subplot(2,1,1);
            plot(timeSpan,dspNodal(node*dpn-dpn+dofToPlot,:));
            title('Displacement (Selected DOF)');
            % Печать спектра.
            RealFFT = fft(dspNodal(node*dpn-dpn+dofToPlot,:));
            % Число семплов.
            nsmp = length(RealFFT);
            % Частота семплирования.
            smplFreq = 1/tStep;
            % Область графика.
            domainFFT = (0:nsmp-1)*smplFreq/nsmp;
            % Модуль.
            absFFT = abs(RealFFT);
            subplot(2,1,2);
            plot(domainFFT,absFFT);
            title('Spectrum');
            hold off;
        end
    end
end

