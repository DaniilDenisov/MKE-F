function test_mises_exports()
% Existing solvers/exporters supply all inputs required by browser recovery.
% Generate real static fixtures for all five forms, including admitted thin walls.
types={'circle','circular_tube','rectangle','square_tube','rectangular_tube'};
base3=fileread(exampleCasePath('CaseSpatialFrame.txt'));
base2=fileread(exampleCasePath('CaseSectionBeam.txt'));
if ~exist('output','dir'), mkdir('output'); end
for i=1:numel(types)
    section=struct('id',1,'name',types{i},'type',types{i}, ...
        'dimensions',struct('b',.1,'h',.2,'d',.1,'t',.003));
    texts={regexprep(base2,'(?m)^\{[^\n]*\}',jsonencode(section)), ...
        [base3 sprintf('\nsections\n1\n%s\nelement_sections\n1\n1,1,27\n',jsonencode(section))]};
    for dimension=2:3
        name=sprintf('Mises%dD_%s',dimension,types{i});
        file=fullfile('output',[name '.txt']); fid=fopen(file,'w');fprintf(fid,'%s',texts{dimension-1});fclose(fid);
        problem=StructFEProblem(file);model=problem.GetAnalysisModel();result=problem.RunSelected();
        data=exportPostprocessorData(model,result,fullfile('output',[name '.json']));
        assert(data.version==dimension+4);
        assert(strcmp(data.model.sections(1).type,types{i}));
        assert(norm(result.equilibriumResidual,inf)<1e-7);
        if dimension==3, assert(isequal(size(result.elementResults(1).sectionForceCoefficients),[6 4])); end
    end
end
fprintf('Static Mises input fixtures: five shapes, existing 2D/v6 and 3D/v7 contracts PASS\n');
end
