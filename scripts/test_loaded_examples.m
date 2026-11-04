function test_loaded_examples()
% Solve the actual TXT downloads produced by test_examples_browser.cjs.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup();
folder = fullfile(root, 'output', 'load-examples');
geometries = {'2d-cantilever', '2d-triangle', '3d-cantilever', '3d-tripod'};
analyses = {'static', 'modal', 'transient'};
for a = 1:numel(analyses)
    for g = 1:numel(geometries)
        name = [geometries{g} '-' analyses{a}];
        problem = StructFEProblem(fullfile(folder, [name '.txt']), struct('verbose',false,'plotting',false));
        model = problem.GetAnalysisModel();
        result = problem.RunSelected();
        assert(strcmp(result.analysisType, analyses{a}));
        if strcmp(analyses{a}, 'modal')
            assert(all(isfinite(result.frequenciesHz)) && all(result.frequenciesHz > 0));
            assert(all(isfinite(result.modeShapes(:))));
        else
            assert(all(isfinite(result.displacements(:))));
            assert(norm(result.displacements, inf) > 0);
            if strcmp(analyses{a}, 'transient')
                % The solver floors duration / timeStep; floating-point division
                % can leave the last sample one step before the requested end.
                assert(numel(result.time) >= 500 && numel(result.time) <= 501);
                assert(all(diff(result.time) > 0));
            else
                assert(norm(result.equilibriumResidual, inf) < 1e-7);
                if ~isempty(strfind(name, 'cantilever'))
                    % Independent Euler-Bernoulli tip deflection, P L^3 / (3 E I).
                    expected = -1000 * 2^3 / (3 * 2e11 * (.1 * .2^3 / 12));
                    assert(abs(result.displacements(model.dofMap(3,2)) - expected) < 1e-12);
                end
            end
        end
        exportPostprocessorData(model, result, fullfile(folder, [name '.json']));
        fprintf('PASS: %s solved and exported\n', name);
    end
end
end
