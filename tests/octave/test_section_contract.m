function test_section_contract()
% Solver provenance, orientations, manual torsion and fresh cross-language cases.
destination=fullfile('output','section-contract');if ~exist(destination,'dir'),mkdir(destination);end
types={'rectangle','circle','square_tube','rectangular_tube','circular_tube','channel','tee','ibeam'};
base2=fileread(exampleCasePath('CaseSectionBeam.txt'));base3=fileread(exampleCasePath('CaseSpatialFrame.txt'));
for dimension=2:3
    for k=1:numel(types)
        if dimension==3 && k>5,continue;end
        s=struct('id',1,'name',types{k},'type',types{k},'dimensions',struct('b',.1,'h',.2,'d',.1,'t',.005,'tw',.008,'tf',.012));
        if k==3 || k==4,s.dimensions.inner_radius=.005;end
        angles=[0 90 180 270];if dimension==3,angles=[0 -37.5];end
        for angle=angles
            if dimension==2
                source=regexprep(base2,'(?m)^\{[^\n]*\}',jsonencode(s));
                source=regexprep(source,'(?m)^1,1,0\r?$',sprintf('1,1,%g',angle));
            else
                if angle~=0,s.torsion=struct('mode','manual','value',.000123);end
                source=[base3 sprintf('\nsections\n1\n%s\nelement_sections\n1\n1,1,%g\n',jsonencode(s),angle)];
            end
            name=sprintf('%dD-%s-%g',dimension,types{k},angle);file=fullfile(destination,[name '.txt']);
            fid=fopen(file,'w');fprintf(fid,'%s',source);fclose(fid);
            problem=StructFEProblem(file,struct('verbose',false,'plotting',false));m=problem.GetAnalysisModel();r=problem.RunSelected();
            data=exportPostprocessorData(m,r,fullfile(destination,[name '.json']));
            assert(data.model.elements(1).sectionAngle==angle);
            p=data.model.elements(1).sectionProperties;sourceProperties=m.elementData(1).sectionProperties;
            for key=fieldnames(p).',assert(isequal(p.(key{1}),sourceProperties.(key{1})));end
            assert(isequal(data.analysis.elementResults(1).sectionStress,r.elementResults(1).sectionStress));
            assert(~isfield(p,'rectangles') && ~isfield(p,'roundedTube'));
            if dimension==3 && angle~=0,assert(p.torsionConstant==.000123 && strcmp(p.torsionMethod,'manual'));end
        end
    end
end
fprintf('  Section contract: 42 fresh solver cases, accepted properties and extrema PASS\n');
end
