function export_mpc_examples()
%EXPORT_MPC_EXAMPLES Generate browser integration artifacts, without editing examples.
root = fileparts(fileparts(mfilename('fullpath')));
setup();
destination = fullfile(root,'output');
if ~exist(destination,'dir'), mkdir(destination); end
cases = {'CaseMPCFrame','CaseMPCModal','CaseMPCTransient','CaseMPCUniform','CaseMPCTruss'};
for i=1:numel(cases)
    p=StructFEProblem(fullfile(root,'examples','cases',[cases{i} '.txt']));
    result=p.RunSelected();
    exportPostprocessorData(p.GetAnalysisModel(),result,fullfile(destination,[cases{i} '.json']));
    if strcmp(result.analysisType,'transient')
        exportPostprocessorData(p.GetAnalysisModel(),result,fullfile(destination,'CaseMPCPartial.json'), ...
            struct('selectedGlobalDOFs',[2 8],'timeStride',3));
    end
end
end
