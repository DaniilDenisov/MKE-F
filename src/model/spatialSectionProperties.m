function p = spatialSectionProperties(s)
% Principal centroidal Y/Z properties. J is stiffness, never mass inertia.
[p,p.rectangles]=sectionProperties(s); d=s.dimensions; p.approximate=false;
mode='auto';
if isfield(s,'torsion')
    if ~isstruct(s.torsion) || ~isscalar(s.torsion) || ~isfield(s.torsion,'mode') || ~ismember(s.torsion.mode,{'auto','manual'})
        error('MKEF:InvalidSection','Invalid torsion mode.');
    end
    mode=s.torsion.mode;
end
if strcmp(mode,'manual')
    if ~isfield(s.torsion,'value'), error('MKEF:InvalidSection','Manual J is required.'); end
    J=s.torsion.value; p.torsionMethod='manual';
elseif ismember(s.type,{'circle','circular_tube'})
    J=p.Iy+p.Iz; p.torsionMethod='exact circular';
elseif strcmp(s.type,'rectangle')
    % Saint-Venant odd series: doi:10.1186/s10033-018-0214-9.
    % The omitted odd tail is bounded by the larger all-integer tail.
    a=max(d.b,d.h); b=min(d.b,d.h); sumOdd=0; n=1; factor=192*b/(pi^5*a);
    while true
        sumOdd=sumOdd+tanh(n*pi*a/(2*b))/n^5;
        bracket=1-factor*sumOdd;
        if factor/(4*n^4)<=1e-12*bracket, break; end
        n=n+2;
    end
    J=a*b^3/3*bracket; p.torsionMethod='Saint-Venant series (relative tail <= 1e-12)';
else
    % https://docs.oasys-software.com/structural/gsa/references-theory/other-sections/
    p.approximate=true;
    if ismember(s.type,{'square_tube','rectangular_tube'})
        h=d.b; if strcmp(s.type,'rectangular_tube'), h=d.h; end
        rm=0; if isfield(p,'roundedTube'), rm=d.inner_radius+d.t/2; end
        enclosed=(d.b-d.t)*(h-d.t)-(4-pi)*rm^2;
        perimeter=2*(d.b+h-2*d.t)-(8-2*pi)*rm;
        J=4*enclosed^2*d.t/perimeter; p.torsionMethod='thin-wall closed';
        if rm>0, p.torsionMethod=[p.torsionMethod ' · rounded median line']; end
    else
        nf=2; if strcmp(s.type,'tee'), nf=1; end
        J=(nf*d.b*d.tf^3+(d.h-nf*d.tf)*d.tw^3)/3; p.torsionMethod='thin-wall open';
    end
end
if ~isnumeric(J) || ~isreal(J) || ~isscalar(J) || ~isfinite(J) || J<=0
    error('MKEF:InvalidSection','J must be positive and finite.');
end
p.torsionConstant=J; p.open=ismember(s.type,{'channel','tee','ibeam'});
end
