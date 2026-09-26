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
        % Вывод диагностической информации в командное окно.
        verbose=true
        % Построение сетки и графиков результатов.
        plotting=true
    end
    methods
        % Конструктор с аргументом.
        function obj = StructFEProblem(filename, options)
            if nargin < 2
                options = struct();
            end
            if ~isstruct(options)
                error('StructFEProblem:InvalidOptions', ...
                    'Options must be provided as a struct.');
            end
            if isfield(options, 'verbose')
                obj.verbose = logical(options.verbose);
            end
            if isfield(options, 'plotting')
                obj.plotting = logical(options.plotting);
            end
            obj.filename = filename;
            % Чтение кейса, создание сетки.
            obj.mesh = FEMesh(obj.filename);
            % Вывод сетки на график.
            if obj.plotting
                obj.mesh.Plot2DMesh();
            end
            % Печать номеров узлов в элементах.
            if obj.verbose
                obj.mesh.DispNN();
            end
            % Печать матрицы соответствия.
            if obj.verbose
                obj.mesh.DispIM();
            end
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
        % Совместимый метод возвращает разбиение СС, не изменяя МЖ и ММ.
        function [fixedDOFs, freeDOFs] = ApplyFixBC(this)
            model = this.GetAnalysisModel();
            [fixedDOFs, freeDOFs] = partitionDOFs(model);
        end
        % Совместимый метод собирает свежий вектор или историю нагрузок.
        function ApplyForceBC(this)
            model = this.GetAnalysisModel();
            if this.ts == 0
                this.F = buildStaticLoad(model);
            else
                if this.tsNum < 1
                    error('MKEF:InvalidTransientOptions', ...
                        'A positive tsNum is required to build a load history.');
                end
                this.F = buildTransientLoad(model, this.ts, this.tsNum);
            end
        end
        % Возвращает копию данных модели для чистого численного ядра.
        function model = GetAnalysisModel(this)
            model = createAnalysisModel(this.K, this.M, this.mesh);
        end
        % Метод запуска расчета статического нагружения.
        function result = RunStatic(this)
            model = this.GetAnalysisModel();
            result = solveStatic(model);
            if this.verbose
                disp('DOFs (displ. components):')
                disp(result.displacements);
                disp('Reactions (components):')
                disp(result.reactions);
            end
        end
        % Метод запуска расчета собственных колебаний.
        function result = RunModal(this)
            model = this.GetAnalysisModel();
            result = solveModal(model);
            if this.verbose
                disp('Natural frequencies (Hz):')
                disp(result.frequenciesHz);
            end
        end
        % Метод запуска анализа динамики со внешними силами.
        function result = RunTransient(this,tStep,tDur,node,dofToPlot)
            options = struct('timeStep', tStep, 'duration', tDur);
            model = this.GetAnalysisModel();
            result = solveTransient(model, options);
            if this.plotting
                globalDOF = this.mesh.iMnod(node, dofToPlot);
                plotTransientResult(result, globalDOF);
            end
        end
    end
end
