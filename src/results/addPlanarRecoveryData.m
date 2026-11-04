function data = addPlanarRecoveryData(data, model, result, settings)
%ADDPLANARRECOVERYDATA V8 planar contract on the existing analysis payload.
data.version=8;
data.metadata.exporter=struct('name','mkef-octave','contractVersion',1);
d=model.dofPerNode; n=2*d;
formulation='linear-bar';if d==3,formulation='euler-bernoulli';end
data.model.numberOfDOFs=model.numberOfDOFs;
data.model.recovery=struct('version',1,'coordinates','local-physical-end', ...
    'stiffnessFormulation',formulation,'massFormulation','consistent-v1');
for i=1:numel(model.elementData)
    e=model.elementData(i);
    check(e.transformation,[n n],'localFromSystem');
    if numel(e.dofs)~=n || any(e.dofs<1 | e.dofs>model.numberOfDOFs | e.dofs~=fix(e.dofs)) || numel(unique(e.dofs))~=n
        error('MKEF:PostprocessorInvalidModel','Invalid planar element DOF mapping.');
    end
    if ~isscalar(e.length) || ~isfinite(e.length) || e.length<=0
        error('MKEF:PostprocessorInvalidModel','Invalid planar element length.');
    end
    data.model.elements(i).length=e.length;
    data.model.elements(i).globalDOFs=e.dofs(:).';
    data.model.elements(i).localFromSystem=e.transformation;
    if strcmp(result.analysisType,'static')
        r=result.elementResults(i);
        check(r.sectionForceCoefficients,[3 4],'sectionForceCoefficients');
        check(r.localDisplacements(:),[n 1],'localDisplacements');
        check(r.equivalentLocalLoadVector(:),[n 1],'equivalentLocalLoadVector');
        data.analysis.elementResults(i).sectionForceCoefficients=r.sectionForceCoefficients;
        data.analysis.elementResults(i).localDisplacements=r.localDisplacements(:).';
        data.analysis.elementResults(i).equivalentLocalLoadVector=r.equivalentLocalLoadVector(:).';
        data.analysis.elementResults(i).sectionStress=r.sectionStress;
    end
end
if strcmp(result.analysisType,'modal')
    values=result.eigenvalues(:).';
    check(values,[1 numel(data.analysis.frequenciesHz)],'eigenvalues');
    data.analysis.eigenvalues=values;
end
if ~strcmp(result.analysisType,'static')
    data.model=addRecoveryOperators(data.model,model,settings.operatorStorage);
end
end
function check(A,shape,name)
if ~isnumeric(A) || ~isreal(A) || ~isequal(size(A),shape) || any(~isfinite(A(:)))
    error('MKEF:PostprocessorInvalidResult','Invalid planar recovery data: %s.',name);
end
end
