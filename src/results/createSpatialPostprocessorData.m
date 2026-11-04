function data = createSpatialPostprocessorData(m,r,settings,a)
%CREATESPATIALPOSTPROCESSORDATA v8; system coordinates may be rotational bases.
% nodeKinematics is sparse COO, one-based [physicalRow,systemDOF,coefficient].
% dofMap columns 4..6 are active-basis slots, NOT necessarily thetaX/Y/Z.
% element.localFromSystem maps its globalDOFs into physical local end DOFs.
required={'numberOfNodes','numberOfDOFs','dofPerNode','dofLabels','nodeCoordinates','elementData', ...
    'dofMap','dofRegistry','releases','warnings','fixedBoundaryConditions','forceBoundaryConditions','multiPointConstraints'};
if ~all(isfield(m,required)), invalid('Incomplete spatial model.'); end
n=m.numberOfDOFs; nn=m.numberOfNodes; d=m.dofPerNode;
if ~isscalar(n) || n<1 || n~=fix(n) || ~isscalar(nn) || nn<1 || nn~=fix(nn) || ...
        ~ismember(d,[3 6]) || ~isequal(size(m.nodeCoordinates),[nn 3]) || isempty(m.elementData) || ...
        ~isequal(size(m.dofMap),[nn d]) || numel(m.dofRegistry)~=n || ~isequal([m.dofRegistry.id],1:n)
    invalid('Inconsistent spatial model dimensions or registry.');
end
code=m.elementData(1).type; spec=elementTypeInfo(code);
if ~ismember(code,[1133 1436]) || spec.dofPerNode~=d || ~isequal(spec.dofLabels(:),m.dofLabels(:)) || any([m.elementData.type]~=code)
    invalid('Inconsistent spatial element family.');
end
baseMap=reshape(1:nn*d,d,nn).'; rebuilt=m.elementData;
for i=1:numel(rebuilt)
    e=m.elementData(i); p=e.properties;
    if code==1436, p=[p e.referenceY e.sectionAngle]; end
    if any(e.nodeNumbers<1) || any(e.nodeNumbers>nn) || any(e.nodeNumbers~=fix(e.nodeNumbers)), invalid('Invalid node reference.'); end
    rebuilt(i)=createStructuralElement(code,m.nodeCoordinates(e.nodeNumbers,:),e.nodeNumbers,p,baseMap);
end
if code==1436
    [map,registry,rebuilt,G,H]=buildSpatialDOFs(nn,rebuilt,m.releases,m.nodeCoordinates);
    if ~isfield(m,'nodalKinematics') || ~isequal(size(G),size(m.nodalKinematics)) || norm(G-m.nodalKinematics,inf)>1e-10 || ...
            ~isfield(m,'rigidMotionMap') || ~isequal(size(H),size(m.rigidMotionMap)) || norm(H-m.rigidMotionMap,inf)>1e-10
        invalid('Inconsistent spatial kinematic maps.');
    end
else
    [map,registry,rebuilt]=buildDOFRegistry(nn,d,rebuilt,m.releases,m.dofLabels); G=speye(n);
end
if ~isequal(map,m.dofMap) || ~isequal(registry,m.dofRegistry), invalid('Inconsistent spatial DOF registry.'); end
for i=1:numel(rebuilt)
    if ~isequal(rebuilt(i).dofs,m.elementData(i).dofs) || ~isequal(size(rebuilt(i).transformation),size(m.elementData(i).transformation)) || ...
            norm(rebuilt(i).transformation-m.elementData(i).transformation,inf)>1e-10 || ...
            norm(rebuilt(i).localBasis-m.elementData(i).localBasis,inf)>1e-10 || ...
            ~isequal(rebuilt(i).nodeCoordinates,m.elementData(i).nodeCoordinates) || ...
            abs(rebuilt(i).length-m.elementData(i).length)>1e-12*rebuilt(i).length
        invalid('Inconsistent spatial element mapping.');
    end
end
c=buildConstraintTransform(m); getElementLoadData(m);
model=struct('dimension',3,'dofPerNode',d,'dofLabels',{m.dofLabels},'numberOfDOFs',n, ...
    'dofMap',m.dofMap,'dofMapMeaning','active nodal basis slots; use nodeKinematics for physical components', ...
    'dofRegistry',m.dofRegistry,'supports',m.fixedBoundaryConditions,'nodalLoads',m.forceBoundaryConditions, ...
    'elementLoads',m.elementLoads,'releases',stripSource(m.releases),'warnings',{m.warnings}, ...
    'mpcs',stripSource(m.multiPointConstraints),'dependentDOFs',c.dependentDOFs,'independentDOFs',c.independentDOFs);
[row,col,val]=find(G); model.nodeKinematics=struct('rows',nn*d,'columns',n,'entries',[row col val]);
model.nodes=repmat(struct('id',0,'x',0,'y',0,'z',0),nn,1);
for i=1:nn, model.nodes(i)=struct('id',i,'x',m.nodeCoordinates(i,1),'y',m.nodeCoordinates(i,2),'z',m.nodeCoordinates(i,3)); end
model.elements=struct([]);
model.sections=m.sections;
for i=1:numel(m.elementData)
    e=m.elementData(i); properties=struct();
    for j=1:numel(spec.properties), properties.(spec.properties{j})=e.properties(j); end
    item=struct('id',i,'type',code,'nodeIds',e.nodeNumbers,'properties',properties,'length',e.length, ...
        'localBasis',e.localBasis,'referenceY',e.referenceY,'sectionAngle',e.sectionAngle, ...
        'globalDOFs',e.dofs(:).','localFromSystem',e.transformation,'sectionId',e.sectionId, ...
        'sectionProperties',exportSectionProperties(e));
    if e.sectionId>0
        s=m.sections(find([m.sections.id]==e.sectionId,1)); p=spatialSectionProperties(s);
        expected=p.area; actual=e.properties(1);
        if code==1436, expected=[expected p.Iy p.Iz p.torsionConstant]; actual=[actual e.properties(5:7)]; end
        if any(abs(actual-expected)>1e-10*abs(expected)), invalid('Section properties disagree with catalog.'); end
    end
    model.elements(i)=item;
end
columns=1; selected=1:n; indices=1;
if strcmp(a.type,'static')
    fields={'displacements','loadVector','reactions'};
    for j=1:numel(fields), f=fields{j}; shape(r.(f),[n 1],f); a.(f)=r.(f); end
    shape(r.equilibriumResidual,[6 1],'equilibriumResidual'); a.equilibriumResidual=r.equilibriumResidual;
    if numel(r.elementResults)~=numel(m.elementData), invalidResult('Wrong element result count.'); end
    a.elementResults=struct([]);
    for i=1:numel(r.elementResults)
        e=r.elementResults(i);
        if e.type~=code || e.elementNumber~=i, invalidResult('Inconsistent element result type or ID.'); end
        for f={'localEndForces','localDisplacements','equivalentLocalLoadVector'}, shape(e.(f{1}),[2*d 1],f{1}); end
        item=struct('elementId',i,'type',code,'localEndForces',e.localEndForces,'localDisplacements',e.localDisplacements, ...
            'equivalentLocalLoadVector',e.equivalentLocalLoadVector,'axialStrain',e.axialStrain,'axialStress',e.axialStress,'axialForce',e.axialForce);
        for f={'axialStrain','axialStress','axialForce'}, shape(e.(f{1}),[1 1],f{1}); end
        if code==1436, shape(e.sectionForceCoefficients,[6 4],'sectionForceCoefficients'); item.sectionForceCoefficients=e.sectionForceCoefficients; end
        item.sectionStress=e.sectionStress;
        a.elementResults(i)=item;
    end
elseif strcmp(a.type,'modal')
    columns=numel(r.eigenvalues); shape(r.eigenvalues,[columns 1],'eigenvalues'); a.eigenvalues=r.eigenvalues;
    indices=1:columns;
else
    columns=numel(r.time); indices=1:settings.timeStride:columns;
    if indices(end)~=columns, indices(end+1)=columns; end
    selected=a.globalDOFIds;
end
for f={'supportReactions','mpcForces','mpcMultipliers'}
    name=f{1};
    if strcmp(a.type,'transient') && ~ismember(name,settings.transientFields), continue; end
    rows=n; if strcmp(name,'mpcMultipliers'), rows=size(c.C,1); end
    shape(r.(name),[rows columns],name);
    ids=selected; if strcmp(name,'mpcMultipliers'), ids=1:rows; end
    a.(name)=r.(name)(ids,indices);
end
if ~strcmp(a.type,'transient') && isfield(r,'nodalSupportReactions')
    for f={'nodalSupportReactions','nodalMPCForces'}
        shape(r.(f{1}),[nn*d columns],f{1}); a.(f{1})=r.(f{1});
    end
end
if strcmp(a.type,'transient')
    a.timeStep=r.timeStep;
    available=false(n,1);available(selected)=true;
    physical=find(full(sum(spones(G(:,~available)),2))==0);
    a.physicalDOFIds=physical(:).';
    for pair={{'supportReactions','nodalSupportReactions'},{'mpcForces','nodalMPCForces'}}
        names=pair{1};if ~ismember(names{1},settings.transientFields),continue;end
        if isfield(r,names{2}),values=r.(names{2});else,values=r.(names{1});end
        shape(values,[nn*d columns],names{2});a.(names{2})=values(physical,indices);
    end
    if ismember('spectrum',settings.transientFields)
        % Transform signed histories before FFT, never amplitude spectra.
        [~,amplitude]=computeSingleSidedSpectrum(full(G(physical,:)*r.displacements),r.timeStep);
        a.nodalDisplacementAmplitudeSpectrum=amplitude;
    end
end
data=struct('format','mkef-postprocessor','version',8,'model',model,'analysis',a);
data.model.recovery=struct('version',1,'coordinates','local-physical-end', ...
    'stiffnessFormulation',spec.formulation,'massFormulation','consistent-v1');
% Static forces are already exported. Dynamic recovery uses the exact
% uncondensed local operators used by the solver, once per element.
if ~strcmp(a.type,'static')
    data.model=addRecoveryOperators(data.model,m,settings.operatorStorage);
end
finiteTree(data);
end
function s=stripSource(s)
if isfield(s,'sourceLine'), s=rmfield(s,'sourceLine'); end
end
function shape(v,sz,name)
if ~isnumeric(v) || ~isreal(v) || issparse(v) || ~isequal(size(v),sz) || any(~isfinite(v(:)))
    invalidResult(['Invalid result shape or values: ' name]);
end
end
function finiteTree(v)
if isnumeric(v) && (~isreal(v) || any(~isfinite(v(:))) || issparse(v)), invalid('Non-finite or sparse exported data.');
elseif isstruct(v)
    fields=fieldnames(v); for i=1:numel(v), for j=1:numel(fields), finiteTree(v(i).(fields{j})); end; end
elseif iscell(v), for i=1:numel(v), finiteTree(v{i}); end
end
end
function invalid(message)
error('MKEF:PostprocessorInvalidModel','%s',message);
end
function invalidResult(message)
error('MKEF:PostprocessorInvalidResult','%s',message);
end
