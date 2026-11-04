function [map,registry,elements,G,H] = buildSpatialDOFs(nodeCount,elements,releases,coords)
%BUILDSPATIALDOFS Orthonormal active nodal rotations and private end rotations.
% G maps system coordinates to physical node-major [u,theta]. H maps a rigid
% translation/rotation about the global origin to system coordinates.
count=numel(elements); released=false(count,2,3); names={'Mx','My','Mz'};
for i=1:numel(releases)
    r=releases(i); axis=find(strcmp(names,r.component));
    if ~isscalar(r.elementId) || ~isfinite(r.elementId) || r.elementId~=fix(r.elementId) || r.elementId<1 || r.elementId>count || ...
            ~isscalar(r.end) || ~ismember(r.end,[1 2]) || isempty(axis)
        error('MKEF:InvalidRelease','Expected existing frame,end (1 or 2),Mx/My/Mz.');
    end
    if released(r.elementId,r.end,axis), error('MKEF:InvalidRelease','Duplicate spatial release.'); end
    released(r.elementId,r.end,axis)=true;
end
directions=cell(nodeCount,1); bases=cell(nodeCount,1);
for i=1:count
    for j=1:2
        node=elements(i).nodeNumbers(j); keep=~reshape(released(i,j,:),3,1);
        directions{node}=[directions{node};elements(i).localBasis(keep,:)];
    end
end
map=zeros(nodeCount,6);
registry=struct('id',{},'kind',{},'nodeId',{},'component',{},'elementId',{},'end',{},'direction',{});
labels={'ux','uy','uz','thetaX','thetaY','thetaZ'};
for node=1:nodeCount
    A=directions{node}; Q=zeros(3,0);
    if ~isempty(A)
        [~,S,V]=svd(A,0); k=min(size(S)); singular=diag(S(1:k,1:k)); rankA=sum(singular>1e-10*max(singular));
        if rankA==3, Q=eye(3); else, Q=V(:,1:rankA); end
        for j=1:size(Q,2)
            [~,k]=max(abs(Q(:,j))); if Q(k,j)<0, Q(:,j)=-Q(:,j); end
        end
    end
    bases{node}=Q;
    for j=1:3+size(Q,2)
        id=numel(registry)+1; map(node,j)=id; v=zeros(1,3);
        if j<=3, v(j)=1; name=labels{j};
        else, v=Q(:,j-3).'; name=labels{j}; if size(Q,2)<3, name=sprintf('rotationBasis%d',j-3); end; end
        registry(id)=struct('id',id,'kind','node','nodeId',node,'component',name,'elementId',0,'end',0,'direction',v);
    end
end
private=zeros(count,2,3);
for i=1:count
    for j=1:2
        for k=1:3
            if ~released(i,j,k), continue; end
            id=numel(registry)+1; private(i,j,k)=id;
            registry(id)=struct('id',id,'kind','elementEnd','nodeId',elements(i).nodeNumbers(j), ...
                'component',labels{k+3},'elementId',i,'end',j,'direction',elements(i).localBasis(k,:));
        end
    end
end
n=numel(registry); G=sparse(6*nodeCount,n); H=sparse(n,6);
for node=1:nodeCount
    rows=(node-1)*6+(1:6); ids=map(node,1:3); Q=bases{node};
    G(rows(1:3),ids)=eye(3); H(ids,1:3)=eye(3);
    p=coords(node,:); H(ids,4:6)=[0 p(3) -p(2);-p(3) 0 p(1);p(2) -p(1) 0];
    ids=map(node,4:3+size(Q,2)); G(rows(4:6),ids)=Q; H(ids,4:6)=Q.';
end
for i=1:count
    B=sparse(12,n); R=elements(i).localBasis;
    for j=1:2
        node=elements(i).nodeNumbers(j); nr=(node-1)*6; er=(j-1)*6;
        B(er+(1:3),:)=R*G(nr+(1:3),:);
        for k=1:3
            id=private(i,j,k);
            if id>0, B(er+3+k,id)=1; H(id,4:6)=R(k,:);
            else, B(er+3+k,:)=R(k,:)*G(nr+(4:6),:); end
        end
    end
    ids=find(any(B,1)); T=full(B(:,ids));
    elements(i).dofs=ids(:); elements(i).transformation=T;
    elements(i).stiffness=T.'*elements(i).localStiffness*T;
    elements(i).mass=T.'*elements(i).localMass*T;
end
end
