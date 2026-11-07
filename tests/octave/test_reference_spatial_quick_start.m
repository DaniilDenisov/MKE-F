function test_reference_spatial_quick_start()
% Keep the downloadable guide's values tied to an independent cantilever solution.
problem=StructFEProblem(exampleCasePath('CaseSpatialQuickStart.txt'));
model=problem.GetAnalysisModel(); result=problem.RunSelected();
e=model.elementData(1); F=-1000; L=2; E=200e9; b=.1; h=.2;
I=b*h^3/12; tip=F*L^3/(3*E*I); stress=abs(F*L)*h/(2*I);
assert(e.type==1436 && e.sectionId==1 && e.sectionAngle==0);
assert(abs(result.displacements(model.dofMap(2,2))-tip)<1e-12);
assert(abs(tip-(-.0002))<1e-12);
assert(abs(result.nodalSupportReactions(2)-1000)<1e-8);
assert(abs(result.elementResults.sectionStress.minimum.value+stress)<1e-6);
assert(abs(result.elementResults.sectionStress.maximum.value-stress)<1e-6);
assert(abs(stress-3e6)<1e-6);
assert(result.elementResults.sectionStress.minimum.xi==0);
assert(result.elementResults.sectionStress.maximum.xi==0);
assert(norm(result.equilibriumResidual,inf)<1e-8);
exportPostprocessorData(model,result,fullfile('output','CaseSpatialQuickStart.json'));
end
