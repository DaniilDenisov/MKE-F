function model = addRecoveryOperators(model, source, storage)
%ADDRECOVERYOPERATORS Export solver-owned local K/M, shared by 2D and 3D.
n=2*model.dofPerNode;
for i=1:numel(source.elementData)
    e=source.elementData(i);
    for field={'localStiffness','localMass'}
        if ~isfield(e,field{1})
            error('MKEF:PostprocessorInvalidModel','Missing solver operator %s.',field{1});
        end
        A=e.(field{1});
        if ~isnumeric(A) || ~isreal(A) || issparse(A) || ~isequal(size(A),[n n]) || any(~isfinite(A(:)))
            error('MKEF:PostprocessorInvalidResult','Invalid recovery operator %s.',field{1});
        end
        if norm(A-A.',inf)>1e-12*max(abs(A(:)))
            error('MKEF:PostprocessorInvalidModel','Asymmetric recovery operator %s.',field{1});
        end
    end
    model.elements(i).recoveryOperators=struct('K',e.localStiffness,'M',e.localMass);
end
model=poolRecoveryOperators(model,storage);
end
