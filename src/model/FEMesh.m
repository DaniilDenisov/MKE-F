% Класс КЭ сетки.
% Copyright 2017 Daniil S. Denisov
classdef FEMesh < handle
    properties (Access = public)
        numberOfElems = 0;
        numberOfNodes = 0;
        numberOfForceBCs = 0;
        numberOfFixBCs = 0;
        dofPerNode;
        iMnod;
        allMeshElems;
        allNodes;
        allFixBCs;
        allForceBCs;
    end
    properties (Access = private)
        sourceFilename = '';
        sourceLineNumber = 0;
        elementType = 0;
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
            obj.allFixBCs = zeros(0, 5);
            obj.allForceBCs = zeros(0, 6);

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
                    case 'nodes'
                        obj.readNodes(fid, markerLine);
                    case 'elems_112'
                        obj.readElements(fid, 112, markerLine);
                    case 'elems_113'
                        obj.readElements(fid, 113, markerLine);
                    case 'bcfix'
                        obj.readFixedConditions(fid);
                    case 'bcforce_stat'
                        obj.readLoads(fid, 10, 5, 'static load');
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

        function readFixedConditions(this, fid)
            this.requireElements(this.sourceLineNumber);
            conditionCount = this.readBlockCount(fid, 'bcfix', false);
            conditions = zeros(conditionCount, 5);
            for i = 1:conditionCount
                [line, lineNumber] = this.readDataLine(fid, 'fixed condition');
                values = this.parseNumericRecord( ...
                    line, lineNumber, 5, ...
                    'fixed condition [type,node,0,0,0]');
                if values(1) ~= fix(values(1)) || ~ismember(values(1), 1:4)
                    this.fail('MKEF:MalformedInput', lineNumber, ...
                        'Unsupported fixed-condition type %g; expected 1, 2, 3, or 4.', ...
                        values(1));
                end
                this.validateNodeID(values(2), lineNumber, 'fixed condition');
                if any(values(3:5) ~= 0)
                    this.fail('MKEF:UnsupportedConstraintValue', lineNumber, ...
                        'Only homogeneous zero constraints are supported.');
                end
                conditions(i, :) = values;
            end
            this.allFixBCs = [this.allFixBCs; conditions];
            this.numberOfFixBCs = size(this.allFixBCs, 1);
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
