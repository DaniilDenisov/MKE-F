function result = example_08_pinned_portal_uniform(verbose)
%EXAMPLE_08_PINNED_PORTAL_UNIFORM Three-element portal under a beam UDL.
%   Checks reactions, all end forces, displacements and the beam moment
%   diagram against an independent force-method solution. The supplied
%   diagram neglects axial strain; both that limit and the finite-EA
%   correction are reported. Pass false to suppress diagnostic output.

if nargin < 1, verbose = true; end
if ~isscalar(verbose) || ...
        ~(islogical(verbose) || (isnumeric(verbose) && isfinite(verbose))) || ...
        ~ismember(double(verbose), [0 1])
    error('MKEF:InvalidExampleOption', 'verbose must be a scalar logical value.');
end
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repositoryRoot);
setup();
problem = StructFEProblem(fullfile(repositoryRoot, 'examples', 'cases', ...
    'CasePinnedPortalUniform.txt'), struct('verbose', false, 'plotting', false));
model = problem.GetAnalysisModel();
analysis = solveStatic(model);

% Independent analytical data, in N and m. A1/A2 are cross-section areas;
% I1/I2 are the drawing's J1/J2, with the same E in all three members.
l = 6; h = 3; q = 10000; E = 2e11;
A1 = 0.02; A2 = 0.02; I1 = 8e-5; I2 = 1.6e-4;
k = (I2/I1)*(h/l);
V = q*l/2;
H0 = q*l^2/(4*h*(2*k+3));
% Release D in X and use horizontal compatibility (Castigliano):
% H*(2*h^3/(3*E*I1) + h^2*l/(E*I2) + l/(E*A2))
%   = q*l^3*h/(12*E*I2).
H = q*l^2/(4*h*(2*k+3+3*I2/(A2*h^2)));
expectedReactions = [H; V; 0; 0; 0; 0; 0; 0; 0; -H; V; 0];
expectedEndForces = [V, H, V; -H, V, H; 0, H*h, 0; ...
    -V, -H, -V; H, V, -H; -H*h, -H*h, H*h];
% AB and DC both run upwards; BC runs to the right.
uB = H*l/(2*E*A2);
vB = -V*h/(E*A1);
thetaA = -uB/h + H*h^2/(6*E*I1);
thetaB = -uB/h - H*h^2/(3*E*I1);
expectedDisplacements = [0; 0; thetaA; uB; vB; thetaB; ...
    -uB; vB; -thetaB; 0; 0; -thetaA];
assert(isequal(analysis.fixedDOFs(:), [1; 2; 10; 11]));
assert(numel(model.elementData) == 3 && model.numberOfDOFs == 12);
near(analysis.reactions, expectedReactions, 1e-9, 1e-7, 'reactions');
near([analysis.elementResults.localEndForces], expectedEndForces, ...
    1e-9, 1e-7, 'local end forces');
near(analysis.displacements, expectedDisplacements, 1e-9, 1e-12, 'displacements');
near(analysis.equilibriumResidual, zeros(3, 1), 0, 1e-7, 'global equilibrium');

% A single loaded beam element suffices for exact end forces. Recover the
% quadratic section moment using equilibrium including q, not a straight
% line between end moments. Sagging is positive on BC.
x = linspace(0, l, 25);
beamForces = analysis.elementResults(2).localEndForces;
beamMoment = -beamForces(3) + beamForces(2)*x - q*x.^2/2;
expectedMoment = -H*h + q*x.*(l-x)/2;
near(beamMoment, expectedMoment, 1e-9, 1e-7, 'beam moment diagram');
near(beamMoment([1 13 25]), [-H*h, q*l^2/8-H*h, -H*h], ...
    1e-9, 1e-7, 'corner and midspan moments');
near([analysis.elementResults.axialForce], [-V, -H, -V], ...
    1e-9, 1e-7, 'axial forces');

% Increase ONLY the beam area: the horizontal reaction must approach the
% drawing's inextensible result. Keep a modest factor to avoid ill conditioning.
stiffAxialModel = model;
beam = model.elementData(2);
properties = beam.properties; properties(1) = 100*properties(1);
stiffAxialModel.elementData(2) = createStructuralElement(beam.type, ...
    beam.nodeCoordinates, beam.nodeNumbers, properties, model.dofMap);
[stiffAxialModel.stiffness, stiffAxialModel.mass] = ...
    assembleGlobalMatrices(stiffAxialModel.elementData, model.numberOfDOFs);
stiffAxialAnalysis = solveStatic(stiffAxialModel);
H100 = q*l^2/(4*h*(2*k+3+3*I2/(100*A2*h^2)));
near(stiffAxialAnalysis.reactions(1), H100, 1e-9, 1e-7, 'large-EA reaction');
assert(abs(stiffAxialAnalysis.reactions(1)-H0) < abs(H-H0)/90);

if verbose
    fprintf('\nExample 08: pinned portal, three frame elements, uniform beam load\n');
    fprintf('l = %g m, h = %g m, q = %g N/m, k = %g\n', l, h, q, k);
    fprintf('Quantity                   Drawing (EA -> Inf)    Finite EA / FEM\n');
    fprintf('VA = VD [N]                %18.9f %18.9f\n', V, analysis.reactions(2));
    fprintf('HA = -HD [N]               %18.9f %18.9f\n', H0, analysis.reactions(1));
    fprintf('MB = MC (sagging +) [N m]   %18.9f %18.9f\n', -H0*h, beamMoment(1));
    fprintf('Mmid (sagging +) [N m]      %18.9f %18.9f\n', q*l^2/8-H0*h, beamMoment(13));
    fprintf('Horizontal reaction difference from drawing: %.6f %%\n', 100*(H0-H)/H0);
    fprintf('Checks: reactions, end forces, displacements, moments, equilibrium, EA limit PASS.\n');
end

result = struct('model', model, 'analysis', analysis, 'k', k, ...
    'drawingH', H0, 'expectedH', H, 'expectedReactions', expectedReactions, ...
    'expectedEndForces', expectedEndForces, ...
    'expectedDisplacements', expectedDisplacements, ...
    'beamX', x, 'beamMoment', beamMoment, 'expectedMoment', expectedMoment);
end

function near(actual, expected, rtol, atol, label)
if ~isequal(size(actual), size(expected)) || any(~isfinite(actual(:))) || ...
        any(abs(actual(:)-expected(:)) > atol + rtol*abs(expected(:)))
    error('MKEF:ReferenceExampleFailed', 'Portal verification failed: %s.', label);
end
end
