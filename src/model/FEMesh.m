% Класс КЭ сетки.
% Copyright 2017 Daniil S. Denisov
classdef FEMesh < handle
    properties (Access = public)
        numberOfElems = 0;
        numberOfNodes = 0;
        numberOfForceBCs = 0;
        numberOfFixBCs = 0;
        dofPerNode;
        numberOfDOFs;
        dofRegistry;
        releases = struct('elementId',{},'end',{},'component',{},'sourceLine',{});
        warnings = {};
        iMnod;
        allMeshElems;
        allNodes;
        allFixBCs;
        allForceBCs;
        elementLoads = struct('type', {}, 'elementId', {}, 'coordinateSystem', {}, 'qx', {}, 'qy', {}, ...
            'qx1', {}, 'qy1', {}, 'qx2', {}, 'qy2', {});
        multiPointConstraints = struct('depNode',{},'depDOF',{},'rhs',{},'masters',{},'sourceLine',{});
        analysisConfiguration;
    end
    properties (Access = private)
        sourceFilename = '';
        sourceLineNumber = 0;
        elementType = 0;
        analysisSourceLine = 0;
    end
    methods
        function obj = FEMesh(filename)
            if isstring(filename) && isscalar(filename)
                filename = char(filename);
            end
            if ~ischar(filename) || isempty(filename)
                error('MKEF:InvalidInputFilename', ...
                    'The input filename must be a non-empty character vector.');
            end

            obj.sourceFilename = filename;
            obj.allNodes = zeros(0, 3);
            obj.allMeshElems = struct([]);
            obj.allFixBCs = normalizeSupports(zeros(0, 5), 3);
            obj.allForceBCs = zeros(0, 6);
            obj.analysisConfiguration = struct();

            [fid, message] = fopen(filename, 'rt');
            if fid < 0
                error('MKEF:InputFileOpenFailed', ...
                    'Cannot open input file "%s": %s', filename, message);
            end
            fileCleanup = onCleanup(@() fclose(fid));

            while true
                [line, reachedEOF] = obj.readLine(fid);
                if reachedEOF
                    break;
                end
                marker = strtrim(line);
                if isempty(marker) || marker(1) == '#'
                    continue;
                end
                markerLine = obj.sourceLineNumber;
                switch marker
                    case 'analysis'
                        obj.readAnalysis(fid, markerLine);
                    case 'nodes'
                        obj.readNodes(fid, markerLine);
                    case 'elems_112'
                        obj.readElements(fid, 112, markerLine);
                    case 'elems_113'
                        obj.readElements(fid, 113, markerLine);
                    case 'mpc'
                        obj.readMPC(fid);
                    case 'releases'
                        obj.readReleases(fid);
                    case 'bcfix'
                        obj.readFixedConditions(fid);
                    case 'bcforce_stat'
                        obj.readLoads(fid, 10, 5, 'static load');
                    case 'eload_uniform'
                        obj.readElementLoads(fid, 20);
                    case 'eload_linear'
                        obj.readElementLoads(fid, 21);
                    case 'bcforce_harm'
                        obj.readLoads(fid, 11, 6, 'harmonic load');
                    case 'bcforce_pulse'
                        obj.readLoads(fid, 12, 5, 'pulse load');
                    case 'bcforce_step'
                        obj.readLoads(fid, 13, 5, 'step load');
                    otherwise
                        obj.fail('MKEF:UnknownInputMarker', markerLine, ...
                            'Unknown section marker "%s".', marker);
                end
            end

            obj.validateCompleteModel();
            [obj.iMnod,obj.dofRegistry,obj.allMeshElems] = buildDOFRegistry( ...
                obj.numberOfNodes,obj.dofPerNode,obj.allMeshElems,obj.releases);
            obj.numberOfDOFs = numel(obj.dofRegistry);
            for i = 1:numel(obj.allFixBCs)
                s = obj.allFixBCs(i);
                if obj.dofPerNode == 3 && s.fixThetaZ && obj.iMnod(s.node,3) == 0
                    message = sprintf('Node %d: thetaZ restraint ignored because all connected ends release Mz.',s.node);
                    obj.warnings{end+1} = message;
                    warning('MKEF:InactiveRotationRestraint','%s',message);
                end
            end
            for i = 1:size(obj.allForceBCs,1)
                node = obj.allForceBCs(i,2);
                if obj.dofPerNode == 3 && obj.iMnod(node,3) == 0 && obj.allForceBCs(i,5) ~= 0
                    error('MKEF:InactiveRotationLoad','Node %d: cannot apply Mz to absent thetaZ.',node);
                end
            end
            obj.validateAnalysisConfiguration();
            validationModel = struct('fixedBoundaryConditions',obj.allFixBCs, ...
                'dofPerNode',obj.dofPerNode,'numberOfNodes',obj.numberOfNodes, ...
                'numberOfDOFs',obj.numberOfDOFs,'dofMap',obj.iMnod, ...
                'multiPointConstraints',obj.multiPointConstraints);
            buildConstraintTransform(validationModel);
            clear fileCleanup;
        end

        function Plot2DMesh(this)
            labels = cellstr(num2str((1:this.numberOfNodes).'));
            plot(this.allNodes(:, 1), this.allNodes(:, 2), 'bo', ...
                'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 2.0);
            hold on;
            grid on;
            for elementNumber = 1:this.numberOfElems
                coordinates = this.allMeshElems(elementNumber).nodeCoordinates;
                plot(coordinates(:, 1), coordinates(:, 2));
                midpoint = mean(coordinates, 1);
                text(midpoint(1), midpoint(2), num2str(elementNumber));
            end
            text(this.allNodes(:, 1), this.allNodes(:, 2), labels, ...
                'LineWidth', 2, 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'right');
            axis equal;
            hold off;
        end

        function DispIM(this)
            format shortG;
            disp('DOF distribution matrix:');
            disp(this.iMnod);
            format compact;
        end

        function DispNN(this)
            format shortG;
            disp('Nodes:');
            for elementNumber = 1:this.numberOfElems
                disp(this.allMeshElems(elementNumber).nodeNumbers);
            end
            format compact;
        end

        function DispElems(this)
            for elementNumber = 1:this.numberOfElems
                element = this.allMeshElems(elementNumber);
                fprintf('type:%d\n', element.type);
                disp(element.nodeCoordinates);
                disp(element.nodeNumbers);
                disp(element.properties);
            end
        end
    end

    methods (Access = private)
        function readAnalysis(this, fid, markerLine)
            if ~isempty(fieldnames(this.analysisConfiguration))
                this.fail('MKEF:DuplicateInputSection', markerLine, ...
                    'The analysis section may appear only once.');
            end

            [line, lineNumber] = this.readDataLine(fid, ...
                'analysis configuration');
            fields = strsplit(line, ',');
            analysisType = lower(strtrim(fields{1}));
            configuration = struct('type', analysisType);

            if ismember(analysisType, {'static', 'modal'})
                if numel(fields) ~= 1
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Analysis %s requires exactly one field.', ...
                        analysisType);
                end
            elseif strcmp(analysisType, 'transient')
                if numel(fields) ~= 5
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        ['Transient analysis requires five fields: ' ...
                         'transient,timeStep,duration,monitorNode,monitorDOF.']);
                end
                values = zeros(1, 4);
                for i = 1:4
                    values(i) = str2double(strtrim(fields{i + 1}));
                end
                if any(~isfinite(values))
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Transient analysis settings must be finite numbers.');
                end
                if values(1) <= 0 || values(2) < values(1)
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        ['Transient timeStep must be positive and duration ' ...
                         'must cover at least one step.']);
                end
                if values(3) ~= fix(values(3)) || values(3) < 1 || ...
                        values(4) ~= fix(values(4)) || values(4) < 1
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Transient monitor node and DOF must be positive integers.');
                end
                configuration.timeStep = values(1);
                configuration.duration = values(2);
                configuration.monitorNode = values(3);
                configuration.monitorDOF = values(4);
            else
                this.fail('MKEF:MalformedInput', lineNumber, ...
                    ['Unsupported analysis type "%s"; expected static, ' ...
                     'modal, or transient.'], analysisType);
            end

            this.analysisConfiguration = configuration;
            this.analysisSourceLine = markerLine;
        end

        function readNodes(this, fid, markerLine)
            if this.numberOfNodes ~= 0
                this.fail('MKEF:DuplicateInputSection', markerLine, ...
                    'The nodes section may appear only once.');
            end
            nodeCount = this.readBlockCount(fid, 'nodes', true);
            nodes = zeros(nodeCount, 3);
            for nodeNumber = 1:nodeCount
                [line, lineNumber] = this.readDataLine(fid, 'node coordinates');
                nodes(nodeNumber, :) = this.parseNumericRecord( ...
                    line, lineNumber, 3, 'node coordinates [X,Y,Z]');
            end
            this.allNodes = nodes;
            this.numberOfNodes = nodeCount;
        end

        function readElements(this, fid, requestedType, markerLine)
            this.requireNodes(markerLine);
            if this.elementType ~= 0 && this.elementType ~= requestedType
                this.fail('MKEF:UnsupportedMixedMesh', markerLine, ...
                    ['Mixed truss/frame meshes are not supported. The file ' ...
                     'already contains type %d elements and then requests type %d.'], ...
                    this.elementType, requestedType);
            end

            elementCount = this.readBlockCount(fid, ...
                sprintf('elems_%d', requestedType), true);
            if this.elementType == 0
                this.elementType = requestedType;
                this.dofPerNode = 2 + (requestedType == 113);
                this.iMnod = reshape( ...
                    1:(this.numberOfNodes*this.dofPerNode), ...
                    this.dofPerNode, this.numberOfNodes).';
            end

            if requestedType == 112
                fieldCount = 6;
                description = ...
                    'truss element [112,node1,node2,A,E,rho]';
            else
                fieldCount = 7;
                description = ...
                    'frame element [113,node1,node2,A,E,rho,I]';
            end
            newElements = cell(elementCount, 1);
            for localNumber = 1:elementCount
                [line, lineNumber] = this.readDataLine(fid, description);
                values = this.parseNumericRecord( ...
                    line, lineNumber, fieldCount, description);
                if values(1) ~= requestedType
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Element type %g does not match section elems_%d.', ...
                        values(1), requestedType);
                end
                nodeNumbers = values(2:3);
                this.validateNodePair(nodeNumbers, lineNumber);
                properties = values(4:end);
                coordinates = this.allNodes(nodeNumbers, :);
                try
                    newElements{localNumber} = createStructuralElement( ...
                        requestedType, coordinates, nodeNumbers, properties, ...
                        this.iMnod);
                catch exception
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Invalid element: %s', exception.message);
                end
            end

            appendedElements = vertcat(newElements{:});
            if isempty(this.allMeshElems)
                this.allMeshElems = appendedElements;
            else
                this.allMeshElems = [this.allMeshElems; appendedElements];
            end
            this.numberOfElems = numel(this.allMeshElems);
        end

        function readReleases(this, fid)
            this.requireElements(this.sourceLineNumber);
            count = this.readBlockCount(fid,'releases',false);
            for i = 1:count
                [line,lineNumber] = this.readDataLine(fid,'release Mz');
                fields = strsplit(line,',');
                if numel(fields) ~= 3
                    this.fail('MKEF:InvalidRelease',lineNumber,'Expected elementId,end,Mz.');
                end
                this.releases(end+1) = struct('elementId',str2double(strtrim(fields{1})), ...
                    'end',str2double(strtrim(fields{2})),'component',strtrim(fields{3}),'sourceLine',lineNumber);
            end
        end

        function readMPC(this, fid)
            this.requireElements(this.sourceLineNumber);
            count = this.readBlockCount(fid, 'mpc', false);
            for i = 1:count
                [line, lineNumber] = this.readDataLine(fid, 'MPC');
                fields = strsplit(line, ',');
                if numel(fields) < 7
                    this.fail('MKEF:MalformedInput',lineNumber,'MPC requires depNode,depDOF,rhs,masterCount,node,dof,coefficient,...');
                end
                masterCount = str2double(fields{4});
                if ~isfinite(masterCount) || masterCount ~= fix(masterCount) || masterCount < 1 || numel(fields) ~= 4+3*masterCount
                    this.fail('MKEF:MalformedInput',lineNumber,'MPC masterCount does not match its fields.');
                end
                v = this.parseNumericRecord(line,lineNumber,4+3*masterCount,'MPC');
                this.multiPointConstraints(end+1) = struct('depNode',v(1),'depDOF',v(2), ...
                    'rhs',v(3),'masters',reshape(v(5:end),3,[]).','sourceLine',lineNumber);
            end
        end

        function readFixedConditions(this, fid)
            this.requireElements(this.sourceLineNumber);
            conditionCount = this.readBlockCount(fid, 'bcfix', false);
            conditions = zeros(conditionCount, 5);
            for i = 1:conditionCount
                [line, lineNumber] = this.readDataLine(fid, 'fixed condition');
                values = this.parseNumericRecord( ...
                    line, lineNumber, 5, ...
                    'fixed condition [type,node,0,0,0]');
                if values(1) ~= fix(values(1)) || ~ismember(values(1), 1:(6 + (this.dofPerNode == 3)))
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Unsupported fixed-condition type %g; expected a supported type 1..7.', ...
                        values(1));
                end
                this.validateNodeID(values(2), lineNumber, 'fixed condition');
                if any(values(3:5) ~= 0)
                    this.fail('MKEF:UnsupportedConstraintValue', lineNumber, ...
                        'Only homogeneous zero constraints are supported.');
                end
                conditions(i, :) = values;
            end
            this.allFixBCs = [this.allFixBCs; normalizeSupports(conditions, this.dofPerNode)];
            this.numberOfFixBCs = size(this.allFixBCs, 1);
        end

        function readElementLoads(this, fid, loadType)
            this.requireElements(this.sourceLineNumber);
            count = this.readBlockCount(fid, 'element load', false);
            fieldCount = 5 + 2*(loadType == 21);
            for i = 1:count
                [line, lineNumber] = this.readDataLine(fid, 'element load');
                v = this.parseNumericRecord(line, lineNumber, fieldCount, 'Element load');
                if v(1) ~= loadType || v(2) ~= fix(v(2)) || v(2) < 1 || ...
                        v(2) > this.numberOfElems || ~ismember(v(3), [1 2]) || ...
                        all(v(4:end) == 0)
                    this.fail('MKEF:InvalidElementLoad', lineNumber, 'Invalid element load.');
                end
                if this.allMeshElems(v(2)).type ~= 113
                    this.fail('MKEF:UnsupportedElementLoad', lineNumber, 'Distributed loads require frame element 113.');
                end
                load = struct('type', loadType, 'elementId', v(2), 'coordinateSystem', v(3), ...
                    'qx', [], 'qy', [], 'qx1', [], 'qy1', [], 'qx2', [], 'qy2', []);
                if loadType == 20
                    load.qx = v(4); load.qy = v(5);
                else
                    load.qx1 = v(4); load.qy1 = v(5); load.qx2 = v(6); load.qy2 = v(7);
                end
                this.elementLoads(end+1) = load;
            end
        end

        function readLoads(this, fid, expectedType, fieldCount, description)
            this.requireElements(this.sourceLineNumber);
            loadCount = this.readBlockCount(fid, description, false);
            loads = zeros(loadCount, 6);
            for i = 1:loadCount
                [line, lineNumber] = this.readDataLine(fid, description);
                values = this.parseNumericRecord( ...
                    line, lineNumber, fieldCount, description);
                if values(1) ~= expectedType
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Load type %g does not match the section type %d.', ...
                        values(1), expectedType);
                end
                this.validateNodeID(values(2), lineNumber, description);
                if this.dofPerNode == 2 && values(5) ~= 0
                    this.fail('MKEF:UnsupportedLoadComponent', lineNumber, ...
                        'A 2D truss node has no Mz degree of freedom.');
                end
                if expectedType == 11 && values(6) <= 0
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Harmonic frequency must be positive.');
                end
                loads(i, 1:fieldCount) = values;
            end
            this.allForceBCs = [this.allForceBCs; loads];
            this.numberOfForceBCs = size(this.allForceBCs, 1);
        end

        function count = readBlockCount(this, fid, blockName, requirePositive)
            [line, lineNumber] = this.readDataLine(fid, [blockName ' count']);
            count = str2double(strtrim(line));
            minimum = 0;
            if requirePositive
                minimum = 1;
            end
            if ~isscalar(count) || ~isfinite(count) || count ~= fix(count) || ...
                    count < minimum
                this.fail('MKEF:MalformedInput', lineNumber, ...
                    'Section %s requires an integer count >= %d.', ...
                    blockName, minimum);
            end
        end

        function values = parseNumericRecord(this, line, lineNumber, ...
                expectedCount, description)
            fields = strsplit(line, ',');
            if numel(fields) ~= expectedCount
                this.fail('MKEF:MalformedInput', lineNumber, ...
                    '%s requires exactly %d comma-separated fields; found %d.', ...
                    description, expectedCount, numel(fields));
            end
            values = zeros(1, expectedCount);
            for i = 1:expectedCount
                values(i) = str2double(strtrim(fields{i}));
            end
            if any(~isfinite(values))
                this.fail('MKEF:MalformedInput', lineNumber, ...
                    '%s contains a missing, nonnumeric, or non-finite value.', ...
                    description);
            end
        end

        function validateNodePair(this, nodeNumbers, lineNumber)
            for i = 1:2
                this.validateNodeID(nodeNumbers(i), lineNumber, 'element');
            end
            if nodeNumbers(1) == nodeNumbers(2)
                this.fail('MKEF:MalformedInput', lineNumber, ...
                    'An element must connect two distinct node IDs.');
            end
        end

        function validateNodeID(this, nodeNumber, lineNumber, description)
            if nodeNumber ~= fix(nodeNumber) || nodeNumber < 1 || ...
                    nodeNumber > this.numberOfNodes
                this.fail('MKEF:MalformedInput', lineNumber, ...
                    '%s refers to invalid node ID %g; valid IDs are 1..%d.', ...
                    description, nodeNumber, this.numberOfNodes);
            end
        end

        function requireNodes(this, lineNumber)
            if this.numberOfNodes == 0
                this.fail('MKEF:InputSectionOrder', lineNumber, ...
                    'The nodes section must appear before element sections.');
            end
        end

        function requireElements(this, lineNumber)
            if this.numberOfElems == 0
                this.fail('MKEF:InputSectionOrder', lineNumber, ...
                    'An element section must appear before boundary conditions.');
            end
        end

        function validateCompleteModel(this)
            lineNumber = max(this.sourceLineNumber, 1);
            if this.numberOfNodes == 0
                this.fail('MKEF:IncompleteInput', lineNumber, ...
                    'The input file contains no nodes section.');
            end
            if this.numberOfElems == 0
                this.fail('MKEF:IncompleteInput', lineNumber, ...
                    'The input file contains no element section.');
            end
        end

        function validateAnalysisConfiguration(this)
            if isempty(fieldnames(this.analysisConfiguration))
                return;
            end

            configuration = this.analysisConfiguration;
            if ~strcmp(configuration.type, 'static') && ~isempty(this.elementLoads)
                this.fail('MKEF:AnalysisLoadMismatch', this.analysisSourceLine, ...
                    'Element loads are supported only in static analysis.');
            end
            if strcmp(configuration.type, 'transient')
                if configuration.monitorNode > this.numberOfNodes
                    this.fail('MKEF:MalformedInput', ...
                        this.analysisSourceLine, ...
                        'Transient monitor node %d is outside 1..%d.', ...
                        configuration.monitorNode, this.numberOfNodes);
                end
                if configuration.monitorDOF > this.dofPerNode
                    this.fail('MKEF:MalformedInput', ...
                        this.analysisSourceLine, ...
                        'Transient monitor DOF %d is outside 1..%d.', ...
                        configuration.monitorDOF, this.dofPerNode);
                end
                if this.iMnod(configuration.monitorNode,configuration.monitorDOF) == 0
                    this.fail('MKEF:InactiveRotationMonitor',this.analysisSourceLine, ...
                        'Node %d: cannot monitor absent thetaZ.',configuration.monitorNode);
                end
            end

            loadTypes = this.allForceBCs(:, 1);
            if strcmp(configuration.type, 'static') && ...
                    any(loadTypes ~= 10)
                this.fail('MKEF:AnalysisLoadMismatch', ...
                    this.analysisSourceLine, ...
                    'Static analysis accepts only load type 10.');
            elseif strcmp(configuration.type, 'modal') && ...
                    ~isempty(loadTypes)
                this.fail('MKEF:AnalysisLoadMismatch', ...
                    this.analysisSourceLine, ...
                    'Modal analysis does not accept nodal load sections.');
            end
        end

        function [line, lineNumber] = readDataLine(this, fid, description)
            while true
                [line, reachedEOF] = this.readLine(fid);
                if reachedEOF
                    this.fail('MKEF:UnexpectedEndOfFile', ...
                        this.sourceLineNumber + 1, ...
                        'Unexpected end of file while reading %s.', description);
                end
                trimmed = strtrim(line);
                if ~isempty(trimmed) && trimmed(1) ~= '#'
                    line = trimmed;
                    lineNumber = this.sourceLineNumber;
                    return;
                end
            end
        end

        function [line, reachedEOF] = readLine(this, fid)
            line = fgetl(fid);
            reachedEOF = isequal(line, -1);
            if ~reachedEOF
                this.sourceLineNumber = this.sourceLineNumber + 1;
            end
        end

        function fail(this, identifier, lineNumber, message, varargin)
            detail = sprintf(message, varargin{:});
            error(identifier, '%s:%d: %s', ...
                this.sourceFilename, lineNumber, detail);
        end
    end
end
