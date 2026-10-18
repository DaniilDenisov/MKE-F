function p = sectionProperties(section, angle)
%SECTIONPROPERTIES Ideal profiles, centroidal axes y (height), z (width).
if nargin < 2, angle=0; end
if ~isscalar(angle) || ~ismember(angle,[0 90 180 270])
    error('MKEF:InvalidSection','Section angle must be 0, 90, 180 or 270.');
end
if ~isstruct(section) || ~isscalar(section) || ~all(isfield(section,{'type','dimensions'})) || ...
        ~ischar(section.type) || ~ismember(section.type,{'rectangle','circle','square_tube','rectangular_tube','circular_tube','channel','tee','ibeam'})
    error('MKEF:InvalidSection','Unsupported section type.');
end
d=section.dimensions;
if strcmp(section.type,'circle'), fields={'d'};
elseif strcmp(section.type,'rectangle'), fields={'b','h'};
elseif strcmp(section.type,'square_tube'), fields={'b','t'};
elseif strcmp(section.type,'rectangular_tube'), fields={'b','h','t'};
elseif strcmp(section.type,'circular_tube'), fields={'d','t'};
else, fields={'b','h','tw','tf'}; end
for i=1:numel(fields)
    k=fields{i};
    if ~isstruct(d) || ~isfield(d,k) || ~isnumeric(d.(k)) || ~isreal(d.(k)) || ~isscalar(d.(k)) || ~isfinite(d.(k)) || d.(k)<=0
        error('MKEF:InvalidSection','Section dimensions must be positive finite numbers.');
    end
end
rects=zeros(0,4); cy=0; cz=0;
if strcmp(section.type,'circular_tube')
    if 2*d.t>=d.d
        error('MKEF:InvalidSection','Section wall thickness leaves no clear space.');
    end
    inner=d.d-2*d.t;
    A=pi*d.t*(d.d-d.t); Iy=A*(d.d^2+inner^2)/16; Iz=Iy;
elseif strcmp(section.type,'circle')
    A=pi*d.d^2/4; Iy=pi*d.d^4/64; Iz=Iy;
else
    if strcmp(section.type,'rectangle')
        rects=[0 0 d.b d.h];
    elseif ismember(section.type,{'square_tube','rectangular_tube'})
        h=d.b; if strcmp(section.type,'rectangular_tube'), h=d.h; end
        if 2*d.t>=min(d.b,h)
            error('MKEF:InvalidSection','Section wall thickness leaves no clear space.');
        end
        rects=[0 0 d.b d.t; 0 h-d.t d.b d.t; ...
            0 d.t d.t h-2*d.t; d.b-d.t d.t d.t h-2*d.t];
    else
        doubleFlange=~strcmp(section.type,'tee');
        if d.tw>=d.b || (1+doubleFlange)*d.tf>=d.h
            error('MKEF:InvalidSection','Section thickness leaves no clear space.');
        end
        rects=[0 d.h-d.tf d.b d.tf];
        if doubleFlange, rects(end+1,:)=[0 0 d.b d.tf]; end
        z=(d.b-d.tw)/2; if strcmp(section.type,'channel'), z=0; end
        rects(end+1,:)=[z doubleFlange*d.tf d.tw d.h-(1+doubleFlange)*d.tf];
    end
    areas=rects(:,3).*rects(:,4); A=sum(areas);
    cy=sum(areas.*(rects(:,2)+rects(:,4)/2))/A;
    cz=sum(areas.*(rects(:,1)+rects(:,3)/2))/A;
    rects(:,1)=rects(:,1)-cz; rects(:,2)=rects(:,2)-cy;
    Iy=sum(rects(:,4).*rects(:,3).^3/12+areas.*(rects(:,1)+rects(:,3)/2).^2);
    Iz=sum(rects(:,3).*rects(:,4).^3/12+areas.*(rects(:,2)+rects(:,4)/2).^2);
end
c=[1 0 -1 0]; s=[0 1 0 -1]; k=angle/90+1;
if isempty(rects), etaMin=-d.d/2; etaMax=d.d/2;
else
    ys=[rects(:,2);rects(:,2);rects(:,2)+rects(:,4);rects(:,2)+rects(:,4)];
    zs=[rects(:,1);rects(:,1)+rects(:,3);rects(:,1);rects(:,1)+rects(:,3)];
    eta=ys*c(k)-zs*s(k); etaMin=min(eta); etaMax=max(eta);
end
if any(~isfinite([A Iy Iz etaMin etaMax])) || any([A Iy Iz]<=0)
    error('MKEF:InvalidSection','Section dimensions overflow or underflow.');
end
p=struct('area',A,'centroidY',cy,'centroidZ',cz,'Iy',Iy,'Iz',Iz, ...
    'momentOfInertia',Iz*c(k)^2+Iy*s(k)^2,'etaMin',etaMin,'etaMax',etaMax);
end
