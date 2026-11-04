function loadHistory = buildTransientLoad(model, timeStep, stepCount)
%BUILDTRANSIENTLOAD Build a fresh load history for one transient run.
% Type 11 follows F0*sin(2*pi*f*t). Type 12 is a one-step pulse at t=dt,
% and type 13 is a persistent step. Columns correspond to t = 0, dt, ..., stepCount*dt.

rejectElementLoads(model);
loadHistory = zeros(model.numberOfDOFs, stepCount + 1);
time = (0:stepCount) * timeStep;
for i = 1:size(model.forceBoundaryConditions, 1)
    boundaryCondition = model.forceBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    if boundaryType==10
        error('MKEF:AnalysisLoadMismatch','Nodal load type 10 is static only; use type 12 for a one-step pulse.');
    end
    nodeNumber = boundaryCondition(2);
    if isfield(model,'dimension') && model.dimension==3
        if isfield(model,'nodalKinematics')
            f=spatialNodalLoad(model,boundaryCondition);
        else
            f=zeros(model.numberOfDOFs,1);
            f(model.dofMap(nodeNumber,:))=getNodalLoadComponents(model,boundaryCondition);
        end
        profile=zeros(1,stepCount+1);
        if boundaryType==12, profile(2)=1;
        elseif boundaryType==11
            frequency=boundaryCondition(9);
            if ~isfinite(frequency) || frequency<=0, error('MKEF:InvalidNodalLoad','Invalid harmonic frequency.'); end
            profile=sin(2*pi*frequency*time);
        elseif boundaryType==13, profile(:)=1;
        else, error('MKEF:UnsupportedLoadType','Unsupported spatial load type.'); end
        loadHistory=loadHistory+f*profile;
        continue;
    end
    nodeDOFs = model.dofMap(nodeNumber, :);
    components = getNodalLoadComponents(model, boundaryCondition);
    if any(components(nodeDOFs == 0) ~= 0)
        error('MKEF:InactiveRotationLoad','Node %d: cannot apply Mz to absent thetaZ.',nodeNumber);
    end
    components = components(nodeDOFs > 0);
    nodeDOFs = nodeDOFs(nodeDOFs > 0);

    if boundaryType == 12
        loadHistory(nodeDOFs, 2) = ...
            loadHistory(nodeDOFs, 2) + components;
    elseif boundaryType == 11
        frequency = boundaryCondition(6);
        loadHistory(nodeDOFs, :) = loadHistory(nodeDOFs, :) + ...
            components * sin(2 * pi * frequency * time);
    elseif boundaryType == 13
        loadHistory(nodeDOFs, :) = loadHistory(nodeDOFs, :) + ...
            repmat(components, 1, stepCount + 1);
    else
        error('MKEF:UnsupportedLoadType', ...
            'Unsupported nodal load type %g.', boundaryType);
    end
end
end
