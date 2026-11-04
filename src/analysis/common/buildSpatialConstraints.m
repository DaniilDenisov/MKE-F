function c = buildSpatialConstraints(model)
%BUILDSPATIALCONSTRAINTS Project physical SPC/MPC onto active coordinates.
% Sparse row elimination avoids a dense null-space matrix. Original rows are
% retained for work-conjugate reaction recovery. No penalty stiffness.
G=model.nodalKinematics; n=size(G,2); np=size(G,1);
physical=rmfield(model,'nodalKinematics'); physical.numberOfDOFs=np;
physical.dofMap=reshape(1:np,6,model.numberOfNodes).';
base=buildConstraintTransform(physical); % Existing topology and cycle checks.
S=sparse(1:numel(base.fixedDOFs),base.fixedDOFs,1,numel(base.fixedDOFs),np);
warnings={}; keep=true(1,size(S,1));
for i=1:size(S,1)
    if norm(S(i,:)*G)<1e-12
        keep(i)=false; warnings{end+1}=sprintf('Physical DOF %d: restraint ignored because its rotation is inactive.',base.fixedDOFs(i));
    end
end
S=S(keep,:);
if isfield(model,'multiPointConstraints')
    index=0;
    for m=model.multiPointConstraints(:).'
        index=index+1;
        refs=[m.depNode m.depDOF; m.masters(:,1:2)];
        contributions=[1;abs(m.masters(:,3))]; referenceScale=0;
        for j=1:size(refs,1)
            projected=G((refs(j,1)-1)*6+refs(j,2),:);
            if norm(projected)<1e-12
                error('MKEF:InvalidMPC','MPC references an inactive physical rotation.');
            end
            referenceScale=max(referenceScale,contributions(j)*norm(projected,inf));
        end
        % Do not normalize cancellation roundoff to unit size. Compare with
        % individual projected terms, so small but active directions remain valid.
        if norm(base.C(index,:)*G,inf)<1e-10*referenceScale
            error('MKEF:InvalidMPC','MPC is redundant after projection onto active rotations.');
        end
    end
end
A=[S; base.C]*G; supportCount=size(S,1);
E=sparse(0,n); pivots=[]; selected=[];
for i=1:size(A,1)
    row=A(i,:); scale=norm(row,inf);
    if scale>0, row=row/scale; end
    for j=1:numel(pivots)
        coefficient=row(pivots(j));
        if coefficient~=0, row=row-coefficient*E(j,:); end
    end
    [maximum,pivot]=max(abs(row));
    if maximum<1e-10
        if i>supportCount, error('MKEF:InvalidMPC','MPC is redundant after projection onto active rotations.'); end
        continue;
    end
    row(abs(row)<1e-13*maximum)=0;
    E(end+1,:)=row/row(pivot); pivots(end+1)=pivot; selected(end+1)=i;
end
independent=setdiff(1:n,pivots); T=sparse(independent,1:numel(independent),1,n,numel(independent));
for j=numel(pivots):-1:1
    row=E(j,:); row(pivots(j))=0; T(pivots(j),:)=-row*T;
end
supportRows=selected(selected<=supportCount); S=S(supportRows,:);
P=[S;base.C]; A=P*G;
fixed=[];
for i=1:size(S,1)
    ids=find(S(i,:)*G); if numel(ids)==1, fixed(end+1)=ids; end
end
c=struct('T',T,'C',base.C*G,'S',S*G,'physicalC',base.C,'physicalS',S, ...
    'A',A,'pivots',pivots,'fixedDOFs',unique(fixed),'freeDOFs',setdiff(1:n,fixed), ...
    'dependentDOFs',pivots(numel(supportRows)+1:end),'independentDOFs',independent, ...
    'order',1:size(base.C,1),'hasMPC',base.hasMPC,'hasTransform',true,'spatial',true,'warnings',{warnings});
if norm(A*T,inf)>1e-9*max(1,norm(A,inf)*norm(T,inf))
    error('MKEF:InvalidConstraint','Spatial constraint elimination failed.');
end
end
