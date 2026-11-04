function [vectors,intensities,resultant] = getSpatialElementLoadData(model)
%GETSPATIALELEMENTLOADDATA Exact integrals of linear full-span force density.
% Intensities: [qx1 qy1 qz1 qx2 qy2 qz2], in the rolled member frame.
n=numel(model.elementData); vectors=cell(n,1); intensities=zeros(n,6); resultant=zeros(6,1);
for i=1:n, vectors{i}=zeros(size(model.elementData(i).localStiffness,1),1); end
if ~isfield(model,'elementLoads') || isempty(model.elementLoads), return; end
for i=1:numel(model.elementLoads)
    q=model.elementLoads(i);
    if ~all(isfield(q,{'type','elementId','coordinateSystem','ends'})) || ...
            ~isscalar(q.elementId) || ~isfinite(q.elementId) || q.elementId~=fix(q.elementId) || q.elementId<1 || q.elementId>n || ...
            ~isscalar(q.type) || ~ismember(q.type,[20 21]) || ~isscalar(q.coordinateSystem) || ~ismember(q.coordinateSystem,[1 2]) || ...
            ~isnumeric(q.ends) || ~isreal(q.ends) || ~isequal(size(q.ends),[3 2]) || any(~isfinite(q.ends(:))) || ~any(q.ends(:))
        error('MKEF:InvalidElementLoad','Invalid spatial distributed load %d.',i);
    end
    e=model.elementData(q.elementId);
    if e.type~=1436, error('MKEF:UnsupportedElementLoad','Distributed loads require a spatial frame.'); end
    if q.type==20 && any(q.ends(:,1)~=q.ends(:,2)), error('MKEF:InvalidElementLoad','Uniform load endpoints differ.'); end
    R=e.localBasis; ends=q.ends;
    if q.coordinateSystem==2, globalEnds=ends; ends=R*ends; else, globalEnds=R.'*ends; end
    a=ends(:,1); b=ends(:,2); L=e.length; f=zeros(12,1);
    f([1 7])=L/6*[2*a(1)+b(1);a(1)+2*b(1)];
    for j=2:3
        if j==2, ids=[2 6 8 12]; sign=1; else, ids=[3 5 9 11]; sign=-1; end
        f(ids)=[L*(7*a(j)+3*b(j))/20;sign*L^2*(3*a(j)+2*b(j))/60; ...
                L*(3*a(j)+7*b(j))/20;-sign*L^2*(2*a(j)+3*b(j))/60];
    end
    vectors{q.elementId}=vectors{q.elementId}+f;
    intensities(q.elementId,:)=intensities(q.elementId,:)+ends(:).';
    force=L*sum(globalEnds,2)/2; firstMoment=L*(globalEnds(:,1)+2*globalEnds(:,2))/6;
    moment=cross(e.nodeCoordinates(1,:).',force)+cross((e.nodeCoordinates(2,:)-e.nodeCoordinates(1,:)).',firstMoment);
    resultant=resultant+[force;moment];
end
end
