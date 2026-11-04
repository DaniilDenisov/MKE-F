function p = exportSectionProperties(e)
%EXPORTSECTIONPROPERTIES Values accepted by the solver, without display geometry.
p=[];
if ~e.sectionId, return; end
p=e.sectionProperties;
for field={'rectangles','roundedTube'}
    if isfield(p,field{1}), p=rmfield(p,field{1}); end
end
keys={'area','centroidY','centroidZ','Iy','Iz','momentOfInertia','etaMin','etaMax'};
for i=1:numel(keys)
    if ~isfield(p,keys{i}) || ~isnumeric(p.(keys{i})) || ~isscalar(p.(keys{i})) || ~isfinite(p.(keys{i}))
        error('MKEF:PostprocessorInvalidModel','Invalid accepted section properties.');
    end
end
expected=p.area; actual=e.properties(1);
if e.type==113, expected=[expected p.momentOfInertia]; actual=[actual e.properties(4)]; end
if e.type==1436, expected=[expected p.Iy p.Iz p.torsionConstant]; actual=[actual e.properties(5:7)]; end
if any(expected<=0) || p.Iy<=0 || p.Iz<=0 || p.momentOfInertia<=0 || p.etaMin>=p.etaMax || any(abs(actual-expected)>1e-10*abs(expected))
    error('MKEF:PostprocessorInvalidModel','Accepted section properties disagree with element.');
end
end
