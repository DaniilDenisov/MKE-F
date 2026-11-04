function out = prepareSpatialJsonData(data)
%PREPARESPATIALJSONDATA Preserve singleton/empty array dimensions in v8 JSON.
out=data; m=data.model;
if isfield(m,'recoveryOperators')
    out.model.recoveryOperators=cell(1,numel(m.recoveryOperators));
    for i=1:numel(m.recoveryOperators)
        op=m.recoveryOperators(i);
        out.model.recoveryOperators{i}=struct('id',op.id,'K',{matrix(op.K)},'M',{matrix(op.M)});
    end
end
for f={'nodes','elements','supports','elementLoads','releases','mpcs','dofRegistry','sections'}
    name=f{1}; out.model.(name)=arrayfun(@(x) x,m.(name),'UniformOutput',false); out.model.(name)=reshape(out.model.(name),1,[]);
end
out.model.dofMap=matrix(m.dofMap); out.model.nodalLoads=matrix(m.nodalLoads);
out.model.nodeKinematics.entries=matrix(m.nodeKinematics.entries);
out.model.dependentDOFs=vector(m.dependentDOFs); out.model.independentDOFs=vector(m.independentDOFs);
for i=1:numel(m.elements)
    for f={'nodeIds','globalDOFs','referenceY'}, out.model.elements{i}.(f{1})=vector(m.elements(i).(f{1})); end
    for f={'localBasis','localFromSystem'}, out.model.elements{i}.(f{1})=matrix(m.elements(i).(f{1})); end
    if isfield(m.elements(i),'recoveryOperators')
        out.model.elements{i}.recoveryOperators.K=matrix(m.elements(i).recoveryOperators.K);
        out.model.elements{i}.recoveryOperators.M=matrix(m.elements(i).recoveryOperators.M);
    end
end
for i=1:numel(m.elementLoads), out.model.elementLoads{i}.ends=matrix(m.elementLoads(i).ends); end
for i=1:numel(m.mpcs), out.model.mpcs{i}.masters=matrix(m.mpcs(i).masters); end
for i=1:numel(m.dofRegistry)
    if isfield(m.dofRegistry,'direction'), out.model.dofRegistry{i}.direction=vector(m.dofRegistry(i).direction); end
end
a=data.analysis;
for f={'frequenciesHz','angularFrequenciesRadPerSec','eigenvalues','equilibriumResidual','time','globalDOFIds','physicalDOFIds','spectrumFrequencyHz'}
    name=f{1}; if isfield(a,name), out.analysis.(name)=vector(a.(name)); end
end
for f={'displacements','loadVector','reactions','supportReactions','mpcForces','mpcMultipliers','nodalSupportReactions','nodalMPCForces','velocities','accelerations','loadHistory','modeShapes','displacementAmplitudeSpectrum','nodalDisplacementAmplitudeSpectrum'}
    name=f{1}; if ~isfield(a,name), continue; end
    if strcmp(a.type,'static'), out.analysis.(name)=vector(a.(name)); else, out.analysis.(name)=matrix(a.(name)); end
end
if isfield(a,'elementResults')
    out.analysis.elementResults=cell(1,numel(a.elementResults));
    for i=1:numel(a.elementResults)
        e=a.elementResults(i);
        for f={'localEndForces','localDisplacements','equivalentLocalLoadVector'}, e.(f{1})=vector(e.(f{1})); end
        if isfield(e,'sectionForceCoefficients'), e.sectionForceCoefficients=matrix(e.sectionForceCoefficients); end
        if isstruct(e.sectionStress), e.sectionStress.positions=vector(e.sectionStress.positions); end
        out.analysis.elementResults{i}=e;
    end
end
end
function out=vector(v)
out=num2cell(double(v(:).'));
end
function out=matrix(v)
out=cell(size(v,1),1); for i=1:size(v,1), out{i}=vector(v(i,:)); end
end
