function stress = recoverSectionStress(p, type, L, f, q)
%RECOVERSECTIONSTRESS Exact normal-stress extrema in normalized coordinate xi.
N=[-f(1) -q(1)*L -(q(3)-q(1))*L/2 0];
M=zeros(1,4);
if type==113, M=[-f(3) f(2)*L q(2)*L^2/2 (q(4)-q(2))*L^2/6]; end
etas=[p.etaMin p.etaMax]; fibers=zeros(2,4); positions=[0 1 unitIntervalRoots(M,1e-8)];
for i=1:2
    fibers(i,:)=N/p.area-M*etas(i)/p.momentOfInertia;
    positions=[positions unitIntervalRoots(fibers(i,2:4).*[1 2 3],1e-8)];
end
% Snap roots within the deduplication tolerance to exact element endpoints.
% Otherwise a root just below 1 can displace the required endpoint itself.
positions(abs(positions)<=1e-10)=0;
positions(abs(positions-1)<=1e-10)=1;
positions=sort(positions); positions=positions([true diff(positions)>1e-10]);
minimum=struct('value',Inf,'xi',0,'eta',0); maximum=minimum; maximum.value=-Inf;
for xi=positions
    for i=1:2
        value=polyval(fliplr(fibers(i,:)),xi);
        point=struct('value',value,'xi',xi,'eta',etas(i));
        if value<minimum.value, minimum=point; end
        if value>maximum.value, maximum=point; end
    end
end
stress=struct('minimum',minimum,'maximum',maximum,'positions',positions);
end
