function test_linear_loads()
%TEST_LINEAR_LOADS Independent analytical beams and virtual-work integration.
p = StructFEProblem(exampleCasePath('CaseUniformFrame.txt'));
m = p.GetAnalysisModel(); L = 2; EI = 2e6;
uniform = m.elementLoads;
for ends = [0 -1000; -1000 0; -300 -900; 600 -600; -1000 -1000].'
    a = ends(1); b = ends(2);
    t = linearModel(m, 0, a, 0, b);
    r = solveStatic(t);
    near(r.reactions(1:3), [0; -L*(a+b)/2; -L^2*(a+2*b)/6]);
    near(r.displacements(5:6), [L^4*(4*a+11*b)/(120*EI); L^3*(a+3*b)/(24*EI)]);
    near(r.elementResults.localEndForces(4:6), zeros(3,1));
    near(r.equilibriumResidual, zeros(3,1));
    s = t; s.fixedBoundaryConditions = [4 1 0 0 0; 6 2 0 0 0];
    rs = solveStatic(s);
    near(rs.reactions([2 5]), -L*[2*a+b; a+2*b]/6);
    near(rs.displacements([3 6]), L^3*[8*a+7*b; -7*a-8*b]/(360*EI));
    near(rs.elementResults.localEndForces([3 6]), [0;0]);
    s.fixedBoundaryConditions = [1 1 0 0 0; 1 2 0 0 0];
    rs = solveStatic(s);
    near(rs.elementResults.localEndForces([3 6]), L^2*[-3*a-2*b;2*a+3*b]/60);
end
t = linearModel(m,0,-1000,0,-1000);
near(solveStatic(t).displacements, solveStatic(m).displacements);
t = linearModel(m,200,0,800,0); r = solveStatic(t);
near(r.displacements(4), L^2*(200+2*800)/(6*.01*2e11));
t = linearModel(m,0,-300,0,-900); r = solveStatic(t);
t.elementLoads(2) = uniform;
near(solveStatic(t).displacements, r.displacements+solveStatic(m).displacements);

% Rotated/translated frame: separate endpoint rotations and independent work.
coords = [3 4; 4.2 5.6];
e = createStructuralElement(113,coords,[1 2],[.01 2e11 7800 1e-5],[1 2 3;4 5 6]);
t = linearModel(m,70,-90,-30,120); t.elementLoads.coordinateSystem = 1;
t.elementData = e; t.nodeCoordinates = coords; t.stiffness = sparse(e.stiffness); t.mass = sparse(e.mass);
[vectors,~,source] = getElementLoadData(t);
near(calculateStaticEquilibrium(t,buildStaticLoad(t),zeros(6,1)),source);
u = [.1;.2;-.3;.4;-.1;.6]; work = 0;
xis = ([-sqrt(3/5) 0 sqrt(3/5)]+1)/2; weights = [5 8 5]/18;
for j = 1:3
    x = xis(j); q = [70;-90]*(1-x)+[-30;120]*x;
    axial = (1-x)*u(1)+x*u(4);
    transverse = (1-3*x^2+2*x^3)*u(2)+L*(x-2*x^2+x^3)*u(3)+ ...
        (3*x^2-2*x^3)*u(5)+L*(-x^2+x^3)*u(6);
    work = work+L*weights(j)*(q(1)*axial+q(2)*transverse);
end
near(vectors{1}.'*u,work);
globalEnds = e.transformation(1:2,1:2).' * [70 -30;-90 120];
t.elementLoads.coordinateSystem = 2;
t.elementLoads.qx1 = globalEnds(1,1); t.elementLoads.qy1 = globalEnds(2,1);
t.elementLoads.qx2 = globalEnds(1,2); t.elementLoads.qy2 = globalEnds(2,2);
[globalVectors,~,globalSource] = getElementLoadData(t);
near(globalVectors{1},vectors{1}); near(globalSource,source);
near(solveStatic(t).equilibriumResidual,zeros(3,1));

source = strrep(fileread(exampleCasePath('CaseLinearFrame.txt')), sprintf('\r'), '');
record = '21,1,2,0,-300,0,-900';
for invalid = {'21,1,2,0,0,0,0','20,1,2,0,1,0,2','21,1,2,0,NaN,0,1', ...
        '21,1,2,0,Inf,0,1','21,1,2,0,,0,1','21,9,2,0,1,0,2','21,1,3,0,1,0,2','21,1,2,0,1,0'}
    mustFail(@() fromText(strrep(source,record,invalid{1})));
end
for task = {'modal','transient,0.01,0.1,2,1'}
    mustFail(@() fromText(strrep(source,sprintf('analysis\nstatic'),sprintf('analysis\n%s',task{1}))));
end
mustFail(@() solveModal(t)); mustFail(@() solveTransient(t,struct()));
bad=t; bad.elementLoads.qx1=NaN; mustFail(@() buildStaticLoad(bad));
bad=t; bad.elementData.type=112; mustFail(@() buildStaticLoad(bad));
mixed = fromText([source sprintf('\neload_uniform\n1\n20,1,1,2,0\neload_linear\n1\n21,1,1,0,1,0,0\n')]);
assert(numel(mixed.mesh.elementLoads)==3);

% Export fixtures also exercise v5 with optional MPC and private rotations.
for name = {'CaseLinearFrame','CaseTriangleFrame','CaseLinearRelease','CaseLinearMPC'}
    p = StructFEProblem(exampleCasePath([name{1} '.txt']));
    model = p.GetAnalysisModel(); r = p.RunStatic();
    near(r.equilibriumResidual,zeros(3,1));
    if strcmp(name{1},'CaseLinearRelease')
        near(r.elementResults.localEndForces([3 6]),[0;0]);
    end
    data = exportPostprocessorData(model,r,fullfile('output',[name{1} '.json']));
    assert(data.version==5);
    json = fileread(fullfile('output',[name{1} '.json']));
    assert(isempty(strfind(json,'"qx":')) && ~isempty(strfind(json,'"qx1":')));
end
file=[tempname() '.json']; cleanup=onCleanup(@() delete(file));
model=mixed.GetAnalysisModel(); exportPostprocessorData(model,mixed.RunStatic(),file);
json=fileread(file); assert(~isempty(strfind(json,'"qx":')) && ~isempty(strfind(json,'"qx1":')));
end

function t = linearModel(m,x1,y1,x2,y2)
t=m; t.elementLoads.type=21; t.elementLoads.qx=[]; t.elementLoads.qy=[];
t.elementLoads.qx1=x1; t.elementLoads.qy1=y1; t.elementLoads.qx2=x2; t.elementLoads.qy2=y2;
end
function near(a,b)
assert(norm(a(:)-b(:),inf) <= 1e-8*max(1,norm(b(:),inf)));
end
function mustFail(action)
try, action(); catch, return; end
error('Expected rejection.');
end
function p = fromText(text)
file=[tempname() '.txt']; fid=fopen(file,'wt'); fprintf(fid,'%s',text); fclose(fid);
cleanup=onCleanup(@() delete(file)); p=StructFEProblem(file);
end
