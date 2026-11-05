function result = recoverSpatialSectionStress(p,coeff)
% sigma_x=N/A+My*z/Iy-Mz*y/Iz; coefficients in ascending powers of x/L.
a=coeff(1,:)/p.area; b=-coeff(6,:)/p.Iz; c=coeff(5,:)/p.Iy;
positions=[0 1]; vertices=[];
curved=isempty(p.rectangles) || isfield(p,'roundedTube');
if curved
    radius=p.etaMax; centers=[0 0];
    if isfield(p,'roundedTube')
        d=p.roundedTube;radius=d.ro;centers=[1 1;1 -1;-1 1;-1 -1].*repmat([d.h/2-radius d.b/2-radius],4,1);
    end
    D=conv(b,b)+conv(c,c); dD=deriv(D);
    for k=1:size(centers,1)
        da=deriv(a+b*centers(k,1)+c*centers(k,2));
        poly=add(4*conv(conv(da,da),D),-radius^2*conv(dD,dD));
        positions=[positions unitIntervalRoots(poly,1e-7) unitIntervalRoots(da,1e-7)];
    end
    positions=[positions unitIntervalRoots(D,1e-7) unitIntervalRoots(b,1e-7) unitIntervalRoots(c,1e-7)];
else
    for k=1:size(p.rectangles,1)
        r=p.rectangles(k,:); vertices=[vertices;r(2) r(1);r(2) r(1)+r(3);r(2)+r(4) r(1);r(2)+r(4) r(1)+r(3)];
    end
    for k=1:size(vertices,1), positions=[positions unitIntervalRoots(deriv(a+b*vertices(k,1)+c*vertices(k,2)),1e-7)]; end
    % Include changes of the supporting corner for exact envelope sampling.
    for k=1:size(vertices,1)
        for j=k+1:size(vertices,1)
            delta=vertices(k,:)-vertices(j,:);
            positions=[positions unitIntervalRoots(b*delta(1)+c*delta(2),1e-7)];
        end
    end
end
positions(abs(positions)<=1e-10)=0; positions(abs(positions-1)<=1e-10)=1;
positions=sort(positions); positions=positions([true diff(positions)>1e-10]);
minimum=struct('value',inf,'xi',0,'y',0,'z',0); maximum=minimum; maximum.value=-inf;
for x=positions
    av=polyval(fliplr(a),x); bv=polyval(fliplr(b),x); cv=polyval(fliplr(c),x);
    points=vertices;
    if isempty(points)
        v=[bv cv]; if norm(v)==0, v=[1 0]; end
        v=radius*v/norm(v);
        if isfield(p,'roundedTube'),v=v+sign(v).*[d.h/2-radius d.b/2-radius];end
        points=[v;-v];
    end
    for k=1:size(points,1)
        y=points(k,1); z=points(k,2); v=struct('value',av+bv*y+cv*z,'xi',x,'y',y,'z',z);
        if v.value<minimum.value, minimum=v; end
        if v.value>maximum.value, maximum=v; end
    end
end
result=struct('minimum',minimum,'maximum',maximum,'positions',positions);
end
function d=deriv(p), d=p(2:end).*(1:numel(p)-1); if isempty(d), d=0; end; end
function v=add(a,b), v=zeros(1,max(numel(a),numel(b))); v(1:numel(a))=a;v(1:numel(b))=v(1:numel(b))+b; end
