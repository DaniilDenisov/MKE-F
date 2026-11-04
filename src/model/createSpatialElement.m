function e = createSpatialElement(type, coords, nodes, p, map)
%CREATESPATIALELEMENT Registry-dispatched spatial elements; local = R*global.
info=elementTypeInfo(type); d=info.dofPerNode;
if ~isnumeric(coords) || ~isreal(coords) || ~isequal(size(coords),[2 3]) || any(~isfinite(coords(:)))
    error('MKEF:InvalidElementGeometry','Spatial elements require two finite XYZ nodes.');
end
delta=coords(2,:)-coords(1,:); L=norm(delta);
if ~isfinite(L) || L<=0, error('MKEF:InvalidElementGeometry','Invalid spatial element length.'); end
if numel(nodes)~=2 || any(nodes~=fix(nodes)) || any(nodes<1) || any(nodes>size(map,1)) || size(map,2)~=d
    error('MKEF:InvalidElementTopology','Invalid spatial connectivity or DOFs.');
end
np=numel(info.properties); extra=4*(type==1436);
if ~isnumeric(p) || ~isreal(p) || numel(p)~=np+extra || any(~isfinite(p)) || any(p(1:np)<=0)
    error('MKEF:InvalidElementProperties','Spatial properties must be positive finite numbers in registry order.');
end
x=delta/L; [~,k]=min(abs(x)); ref=zeros(1,3); ref(k)=1;
angle=0;
if type==1436, ref=p(8:10); angle=p(11); end
if norm(ref)==0 || ~isfinite(norm(ref)), error('MKEF:InvalidElementOrientation','Invalid reference Y vector.'); end
ref=ref/norm(ref); y=ref-dot(ref,x)*x;
if norm(y)<1e-10, error('MKEF:InvalidElementOrientation','Reference Y is parallel or nearly parallel to the member.'); end
y=y/norm(y); z=cross(x,y);
R=[x;cosd(angle)*y+sind(angle)*z;-sind(angle)*y+cosd(angle)*z];
K=zeros(2*d); K([1 d+1],[1 d+1])=p(1)*p(2)/L*[1 -1;-1 1];
T=blkdiag(R,R);
if type==1436
    T=blkdiag(R,R,R,R);
    K([4 10],[4 10])=p(4)*p(7)/L*[1 -1;-1 1];
    B=[12 6*L -12 6*L;6*L 4*L^2 -6*L 2*L^2;-12 -6*L 12 -6*L;6*L 2*L^2 -6*L 4*L^2];
    K([2 6 8 12],[2 6 8 12])=p(2)*p(6)/L^3*B;
    S=diag([1 -1 1 -1]);
    K([3 5 9 11],[3 5 9 11])=p(2)*p(5)/L^3*S*B*S;
end
M=p(3)*p(1)*L/6*kron([2 1;1 2],eye(3));
if type==1436
    M=zeros(12); M([1 7],[1 7])=p(3)*p(1)*L/6*[2 1;1 2];
    % Polar area inertia governs torsional mass, not the torsion constant J.
    M([4 10],[4 10])=p(3)*(p(5)+p(6))*L/6*[2 1;1 2];
    B=p(3)*p(1)*L/420*[156 22*L 54 -13*L;22*L 4*L^2 13*L -3*L^2;54 13*L 156 -22*L;-13*L -3*L^2 -22*L 4*L^2];
    M([2 6 8 12],[2 6 8 12])=B;
    M([3 5 9 11],[3 5 9 11])=S*B*S;
end
e=struct('type',type,'nodeNumbers',nodes(:).','nodeCoordinates',coords, ...
    'properties',reshape(p(1:np),1,[]),'dofs',reshape(map(nodes,:).',[],1),'length',L, ...
    'transformation',T,'localStiffness',K,'localMass',M,'stiffness',T.'*K*T, ...
    'mass',T.'*M*T,'localBasis',R,'sectionId',0,'sectionAngle',angle,'sectionProperties',struct(), ...
    'referenceY',ref);
if any(~isfinite(e.stiffness(:))) || any(~isfinite(e.mass(:)))
    error('MKEF:InvalidElementProperties','Element stiffness or mass overflow.');
end
end
