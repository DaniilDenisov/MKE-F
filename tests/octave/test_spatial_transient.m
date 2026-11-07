function test_spatial_transient()
p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt')); m=p.GetAnalysisModel(); m.forceBoundaryConditions=zeros(0,9);
modal=solveModal(m); phi=modal.modeShapes(:,1); omega=modal.angularFrequenciesRadPerSec(1);
errors=[];
for steps=[100 200]
    duration=1.3*2*pi/omega; dt=duration/steps;
    r=solveTransient(m,struct('timeStep',dt,'duration',duration,'initialDisplacement',phi));
    expected=phi*cos(omega*r.time(end)); errors(end+1)=norm(r.displacements(:,end)-expected)/norm(phi);
    energy=sum(r.velocities.*(m.mass*r.velocities))/2+sum(r.displacements.*(m.stiffness*r.displacements))/2;
    assert(max(abs(energy-energy(1)))/energy(1)<1e-9);
    assert(norm(r.nodalDisplacements-m.nodalKinematics*r.displacements,inf)==0);
end
assert(errors(2)<.3*errors(1));
% Independent closed-form forced SDOF: the one spatial truss free coordinate.
p=StructFEProblem(exampleCasePath('CaseSpatialTruss.txt')); t=p.GetAnalysisModel(); mode=solveModal(clearLoads(t));
w0=mode.angularFrequenciesRadPerSec(1); k=2e11*.01/2; F=1000; dt=2*pi/w0/1000;
t.forceBoundaryConditions=[13 2 0 0 F 0 0 0 0];
r=solveTransient(t,struct('timeStep',dt,'duration',300*dt));
expected=F/k*(1-cos(w0*r.time)); assert(norm(r.displacements(6,:)-expected,inf)<1e-5*F/k);
w=.7*w0; t.forceBoundaryConditions=[11 2 0 0 F 0 0 0 w/(2*pi)];
r=solveTransient(t,struct('timeStep',dt,'duration',300*dt));
expected=F/k/(1-(w/w0)^2)*(sin(w*r.time)-w/w0*sin(w0*r.time));
assert(norm(r.displacements(6,:)-expected,inf)<1e-5*F/k);
t.forceBoundaryConditions(1)=12; pulse=buildTransientLoad(t,dt,5);
assert(nnz(pulse)==1 && pulse(6,2)==F);
% Frequency is column 9, not a moment component.
m.forceBoundaryConditions=[11 2 0 0 0 400 0 0 7];
loads=buildTransientLoad(m,.001,10); assert(norm(loads(10,:)-400*sin(2*pi*7*(0:10)*.001),inf)<1e-12);
% An end release carries no external moment; elastic and inertia moments cancel.
text=fileread(exampleCasePath('CaseSpatialFrame.txt')); text=regexprep(text,'bcforce_stat_3d[\s\S]*$','');
text=[text sprintf('\nreleases\n1\n1,2,Mz\n')]; p=fromText(text); a=p.GetAnalysisModel(); modes=solveModal(a);
r=solveTransient(a,struct('timeStep',.0001,'duration',.01,'initialDisplacement',modes.modeShapes(:,1)));
private=find(strcmp({a.dofRegistry.kind},'elementEnd'));
assert(norm(r.reactions(private,:),inf)<1e-8*norm(a.stiffness*r.displacements,inf));
assert(norm(a.stiffness(private,:)*r.displacements,inf)>1e-6);
% Projected supports/MPC and initial conditions use the same transform.
a.multiPointConstraints=struct('depNode',2,'depDOF',2,'rhs',0,'masters',[2 3 2]);
modes=solveModal(a); r=solveTransient(a,struct('timeStep',.0001,'duration',.005,'initialDisplacement',modes.modeShapes(:,1)));
assert(norm(r.nodalDisplacements(8,:)-2*r.nodalDisplacements(9,:),inf)<1e-12);
bad=zeros(a.numberOfDOFs,1); bad(1)=1; failed=false;
try, solveTransient(a,struct('timeStep',.001,'duration',.01,'initialDisplacement',bad)); catch err, failed=strcmp(err.identifier,'MKEF:InvalidInitialConditions'); end
assert(failed);
end
function m=clearLoads(m)
m.forceBoundaryConditions=zeros(0,9);
end
function p=fromText(text)
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid); cleanup=onCleanup(@() delete(file)); p=StructFEProblem(file);
end
