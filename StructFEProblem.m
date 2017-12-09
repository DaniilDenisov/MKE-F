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
        % Шаг по времени (def=0).
        ts=0
        % Число шагов по времени (def=0).
        tsNum=0
        % Время задачи (def=0).
        tDur=0
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
            % масс в зависимости от кол-ва СС на узел.
            dofPerNode = obj.mesh.dofPerNode;
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
                    GlobM, obj.mesh.iMnod);
            end
            obj.K = GlobK;
            obj.M = GlobM;
            % Вектор правой части (сил). Преаллокация без ГУ.
            obj.F = zeros(obj.mesh.numberOfNodes*dofPerNode,1);
            
        end
        % Метод наложения ГУ заделки. Действует на ММ и МЖ.
        % Годен для статики, динамики и собств. колебаний.
        function ApplyFixBC(this)
            BCs = this.mesh.allFixBCs;
            BCtotal = this.mesh.numberOfFixBCs;
            for i=1:BCtotal
                % Взять ГУ
                currBC = BCs(i,:);
                % Взять тип ГУ.
                typeBC = currBC(1);
                % Взять номер узла.
                nnumBC = currBC(2);
                % Вычислить номер степени свободы в глоб. МЖ, MM, ВПЧ.
                GLDOFs = this.mesh.iMnod(nnumBC,:);
                % Применить ГУ к МЖ и ММ.
                % Заделка всех СС. Обнулить строки, столбцы и записать 1 на диагональ.
                if (typeBC==1)
                    for dofNum=1:this.mesh.dofPerNode
                        % Матрица масс.
                        this.M(GLDOFs(dofNum),:)=0; % Обнулить строки.
                        this.M(:,GLDOFs(dofNum))=0; % Обнулить столбцы.
                        this.M(GLDOFs(dofNum),GLDOFs(dofNum))=1; % Единица.
                        % Матрица жесткости.
                        this.K(GLDOFs(dofNum),:)=0;
                        this.K(:,GLDOFs(dofNum))=0;
                        % Если расчет статический, то МЖ нельзя оставлять
                        % сингулярной (с нулями без единичек на диагонали)
                        if this.ts==0
                            this.K(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        end
                    end
                end
                % Заделка кроме гор. СС, обнулить все, кроме первой СС.
                if (typeBC==2)
                    for dofNum=1:this.mesh.dofPerNode
                        % Пропускать первую СС.
                        if dofNum==1
                            continue
                        end
                        % Матрица масс.
                        this.M(GLDOFs(dofNum),:)=0;
                        this.M(:,GLDOFs(dofNum))=0;
                        this.M(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        % Матрица жесткости.
                        this.K(GLDOFs(dofNum),:)=0;
                        this.K(:,GLDOFs(dofNum))=0;
                        % Если расчет статический, то МЖ нельзя оставлять
                        % сингулярной (с нулями без единичек на диагонали)
                        if this.ts==0
                            this.K(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        end
                    end
                end
                % Заделка кроме верт. СС, обнулить все, кроме второй СС.
                if (typeBC==3)
                    for dofNum=1:this.mesh.dofPerNode
                        % Пропускать вторую СС.
                        if dofNum==2
                            continue
                        end
                        % Матрица масс.
                        this.M(GLDOFs(dofNum),:)=0;
                        this.M(:,GLDOFs(dofNum))=0;
                        this.M(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        % Матрица жесткости.
                        this.K(GLDOFs(dofNum),:)=0;
                        this.K(:,GLDOFs(dofNum))=0;
                        % Если расчет статический, то МЖ нельзя оставлять
                        % сингулярной (с нулями без единичек на диагонали)
                        if this.ts==0
                            this.K(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        end
                    end
                end
                % Заделка всех кроме угла поворота.
                if (typeBC==4)
                    for dofNum=1:this.mesh.dofPerNode
                        % Пропускать третью СС.
                        if dofNum==3
                            continue
                        end
                        % Матрица масс.
                        this.M(GLDOFs(dofNum),:)=0;
                        this.M(:,GLDOFs(dofNum))=0;
                        this.M(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        % Матрица жесткости.
                        this.K(GLDOFs(dofNum),:)=0;
                        this.K(:,GLDOFs(dofNum))=0;
                        % Если расчет статический, то МЖ нельзя оставлять
                        % сингулярной (с нулями без единичек на диагонали)
                        if this.ts==0
                            this.K(GLDOFs(dofNum),GLDOFs(dofNum))=1;
                        end
                    end
                end
            end
        end
        % Метод наложения ГУ усилий.
        function ApplyForceBC(this)
            BCs = this.mesh.allForceBCs;
            BCtotal = this.mesh.numberOfForceBCs;
            for i=1:BCtotal
                % Взять ГУ
                currBC = BCs(i,:);
                % Взять тип ГУ.
                typeBC = currBC(1);
                % Взять номер узла.
                nnumBC = currBC(2);
                % Вычислить номер степеней свободы в глоб. ВПЧ.
                GLDOFs = this.mesh.iMnod(nnumBC,:);
                % Применить ГУ к ВПЧ.
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
                % Если прикладывается сила гармоническая, создать tsNum
                % столбцов с сохранением неизменных сил.
                if typeBC==11
                    % Проверка не вызвана ли ApplyForceBC без tStep.
                    if this.ts==0
                        error('No timestep during harm. BC application!');
                    end
                    forceValue = zeros(3,1);
                    % Считать Fx,Fy,Fz в forceValue.
                    forceValue(1) = currBC(1,3);
                    forceValue(2) = currBC(1,4);
                    forceValue(3) = currBC(1,5);
                    % Считать частоту в forceValue.
                    freq = currBC(1,6);
                    % Формирование вектора правой части для каждого шага
                    % по времени.
                    for s=1:this.tsNum
                        this.F(GLDOFs(1),s) = this.F(GLDOFs(1),s)+...
                            forceValue(1)*sin((2*pi*freq)*(s*this.ts));
                        this.F(GLDOFs(2),s) = this.F(GLDOFs(2),s)+...
                            forceValue(2)*sin((2*pi*freq)*(s*this.ts));
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
            this.ApplyFixBC();
            this.ApplyForceBC();
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
            this.ApplyFixBC();
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
            % Сохранить в поля объекта задачи: Целое число временных
            % шагов, длительность задачи и шаг.
            this.tsNum = fix(tDur/tStep);
            this.tDur = tDur;
            this.ts = tStep;
            % Массивы для узловых скоростей, ускорений и перемещений.
            nn = this.mesh.numberOfNodes;
            dpn = this.mesh.dofPerNode;
            acsNodal = zeros(nn*dpn,1);
            speNodal = zeros(nn*dpn,1);
            dspNodal = zeros(nn*dpn,1);
            % Установка вектора правой части в ноль для всех шагов.
            this.F = zeros(nn*dpn,this.tsNum);
            % ГУ.
            this.ApplyFixBC();
            this.ApplyForceBC();
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
            this.ApplyFixBC();
            % Формирование эффективной матрицы жесткости.
            Khat = this.K+a0*this.M;
            % Для каждого шага по времени считаем перемещения.
            for i=1:this.tsNum
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
            timeSpan = 0:tStep:this.tsNum*tStep;
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

