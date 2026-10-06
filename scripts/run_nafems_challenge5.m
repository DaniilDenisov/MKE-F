function summary = run_nafems_challenge5(profile)
%RUN_NAFEMS_CHALLENGE5 Reproduce saved cases without changing solver or inputs.
if nargin == 0, profile = 'quick'; end
if ~ischar(profile) || ~ismember(profile, {'quick','full'})
    error('NAFEMS5:Profile','Expected quick or full.');
end
root = fileparts(fileparts(mfilename('fullpath'))); setup();
source = fullfile(root,'examples','cases','nafems-challenge-5');
destination = fullfile(root,'output','nafems-challenge-5',profile);
if exist(destination,'dir') ~= 7, mkdir(destination); end
log = fopen(fullfile(destination,'run.log'),'wt'); assert(log>=0);
cleanup = onCleanup(@() fclose(log));
fprintf(log,'status=RUNNING\nprofile=%s\nOctave=%s\nstarted=%s\n',profile,version(),datestr(now,31));
% Run git from the repository even when Octave was launched elsewhere.
previousDir = pwd; cd(root); restore = onCleanup(@() cd(previousDir));
[gitStatus,revision] = system('git rev-parse HEAD');
if gitStatus ~= 0, revision = 'unavailable'; end
[~,dirty] = system('git status --porcelain');
fprintf(log,'revision=%s\nworking_tree_changes=%s\n',strtrim(revision),strtrim(dirty));
clear restore;
fprintf(log,'mass=consistent Euler-Bernoulli; no separate rho*I contribution\n');
fprintf(log,'old artifacts may remain after a failed rerun; only this log lists completed cases\n');
frequencyFile = fopen(fullfile(destination,'frequencies.csv'),'wt'); assert(frequencyFile>=0);
frequencyCleanup = onCleanup(@() fclose(frequencyFile));
fprintf(frequencyFile,'file,mode,frequency_hz,analytical_hz,relative_error\n');
checksFile = fopen(fullfile(destination,'checks.csv'),'wt'); assert(checksFile>=0);
checksCleanup = onCleanup(@() fclose(checksFile));
fprintf(checksFile,'file,independent_dofs,mass_kg,max_relative_residual\n');
summary = struct('profile',profile,'status','INCOMPLETE','output',destination);
currentCase = ''; frameRuns = {}; analyticalRows = [];
try
    fid = fopen(fullfile(source,'manifest.csv'),'rt'); assert(fid>=0);
    manifestCleanup = onCleanup(@() fclose(fid));
    records = textscan(fid,'%s%s%f%s','Delimiter',',','HeaderLines',1);
    clear manifestCleanup;
    converged = false;
    for row = 1:numel(records{1})
        name=records{1}{row}; kind=records{2}{row}; n=records{3}(row); group=records{4}{row};
        if strcmp(profile,'quick') && ~strcmp(group,'quick'), continue; end
        if strcmp(group,'conditional') && converged, continue; end
        currentCase = name;
        fprintf('NAFEMS5 %s: %s\n',profile,name);
        fprintf(log,'START %s\n',name); fflush(log);
        problem = StructFEProblem(fullfile(source,name),struct('verbose',false,'plotting',false));
        model = problem.GetAnalysisModel(); result = problem.RunSelected();
        % Export first: a failed check still leaves inspectable solver output.
        exportPostprocessorData(model,result,fullfile(destination,strrep(name,'.txt','.json')), ...
            struct('title',name,'lengthUnit','m','forceUnit','N','momentUnit','N*m','timeUnit','s'));
        metrics = nafems_challenge5_check(model,result,kind,n);
        fprintf(checksFile,'%s,%d,%.17g,%.17g\n',name,metrics.dofs,metrics.massKg,metrics.residual);
        fflush(checksFile);
        if strcmp(result.analysisType,'modal')
            frequencies=result.frequenciesHz;
            for k=1:numel(frequencies)
                exact=NaN;
                if strcmp(kind,'axial'), exact=(2*k-1)/4*sqrt(210e9/7800); end
                if strcmp(kind,'bending') && k<=2, exact=k^2*pi/2*sqrt(210e9*(pi*.01^4/4)/(7800*pi*.01^2)); end
                fprintf(frequencyFile,'%s,%d,%.17g,%.17g,%.17g\n',name,k,frequencies(k),exact,frequencies(k)/exact-1);
            end
            fflush(frequencyFile);
            if strcmp(kind,'axial') || strcmp(kind,'bending')
                analyticalRows(end+1,:) = [strcmp(kind,'bending') n metrics.error1 metrics.error2];
            elseif strcmp(kind,'frame')
                last = find(frequencies>1600,1);
                if isempty(last)
                    if n>=25, error('NAFEMS5:FrequencyRange','No mode above 1600 Hz in %s.',name); end
                    last=numel(frequencies); % Coarse meshes have only 6*n modes.
                end
                [points,shapes]=nafems_challenge5_shapes(model,result.modeShapes(:,1:last));
                snapshot=struct('n',n,'frequencies',frequencies(1:last),'points',points,'shapes',shapes);
                frameRuns{end+1}=snapshot;
                writeShapes(destination,snapshot);
                if n>=200
                    [converged,comparison]=compareMeshes(frameRuns{end-1},snapshot);
                    writeComparison(destination,comparison,frameRuns{end-1}.n,n);
                    fprintf(log,'CONVERGENCE %d -> %d: %d; count_below_1600=%d\n',frameRuns{end-1}.n,n,converged,last-1);
                end
            end
        end
        fprintf(log,'PASS %s dofs=%d residual=%.6g\n',name,metrics.dofs,metrics.residual); fflush(log);
    end
    if strcmp(profile,'full') && ~converged
        error('NAFEMS5:NoConvergence','The 200 -> 400 comparison did not meet acceptance criteria.');
    end
    makePlots(destination,analyticalRows,frameRuns);
    latest=frameRuns{end}; count=sum(latest.frequencies<=1600);
    report=fopen(fullfile(destination,'comparison.txt'),'wt'); assert(report>=0);
    reportCleanup=onCleanup(@() fclose(report));
    fprintf(report,'MKE-F frame: n=%d, modes <=1600 Hz: %d; NAFEMS reports 22.\n',latest.n,count);
    fprintf(report,'First frequencies: %.9g, %.9g Hz.\n',latest.frequencies(1:2));
    fprintf(report,'NAFEMS figure 6 beam labels: 366, 496, 1021, 1183 Hz; truss: 451, 1027 Hz.\n');
    for target=[451 1027]
        below=find(latest.frequencies<=target,1,'last'); above=find(latest.frequencies>target,1);
        fprintf(report,'Computed beam modes around %g Hz: %.9g (mode %d), %.9g (mode %d) Hz.\n', ...
            target,latest.frequencies(below),below,latest.frequencies(above),above);
    end
    fprintf(report,'These are diagnostic comparisons, not identical-mass reference tests.\n');
    fprintf(report,'Current consistent mass has no separate rho*I term. This identifies a formulation difference, not a proof that it explains every discrepancy.\n');
    fprintf(report,'No parameter fitting or software changes. Figure 7 and lumped mass are outside scope.\n');
    if strcmp(profile,'quick'), fprintf(report,'Quick profile does not certify convergence of the full band.\n'); end
    clear reportCleanup;
    summary.status='PASS'; summary.elementsPerMember=latest.n; summary.modesBelow1600=count;
    fprintf(log,'status=PASS\nfinished=%s\n',datestr(now,31));
catch exception
    fprintf(log,'status=INCOMPLETE\ncase=%s\nerror=%s\nmessage=%s\n',currentCase,exception.identifier,exception.message);
    fprintf(log,'Do not change solver, mass, formats or tolerances to bypass this failure. Investigate and report first.\n');
    fflush(log);
    rethrow(exception);
end
fprintf('NAFEMS5 %s PASS: %s\n',profile,destination);
end

function [passed,rows] = compareMeshes(previous,current)
old=previous.frequencies; new=current.frequencies;
count=sum(new<=1600); rows=zeros(count,5);
passed = sum(old<=1600)==count;
% Quadrature weights on the same 201 points/member; include member length.
w=[ones(201,1);sqrt(2)*ones(201,1)]; w([1 201 202 402])=w([1 201 202 402])/2;
w=kron(w,ones(2,1));
X=previous.shapes.*sqrt(w); Y=current.shapes.*sqrt(w);
X=X./sqrt(sum(X.^2,1)); Y=Y./sqrt(sum(Y.^2,1)); mac=abs(X.'*Y).^2;
used=false(numel(old),1);
for k=1:count
    candidates=find(~used & abs(old/new(k)-1)<.05);
    if isempty(candidates), passed=false; rows(k,:)=[k NaN Inf 0 0]; continue; end
    [score,pos]=max(mac(candidates,k)); j=candidates(pos); used(j)=true;
    difference=abs(new(k)-old(j))/new(k);
    tolerance=.01; if k<=2, tolerance=.001; end
    rows(k,:)=[k j difference score tolerance];
    passed=passed && difference<tolerance;
    % Near-degenerate modes are compared by their weighted subspaces below.
    cluster=find(abs(new/new(k)-1)<1e-4);
    if numel(cluster)==1
        passed=passed && score>.95;
    else
        oldCluster=find(abs(old/new(k)-1)<.05);
        [U,~]=qr(X(:,oldCluster),0); [V,~]=qr(Y(:,cluster),0);
        passed=passed && numel(oldCluster)==numel(cluster) && min(svd(U.'*V))^2>.95;
    end
end
end

function writeComparison(destination,rows,oldN,newN)
fid=fopen(fullfile(destination,sprintf('convergence_n%03d_n%03d.csv',oldN,newN)),'wt');
assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'current_mode,previous_mode,relative_frequency_change,MAC,frequency_tolerance\n');
fprintf(fid,'%d,%g,%.17g,%.17g,%.17g\n',rows.');
end

function selected = selectedModes(frequencies)
selected=1;
for target=[451 1027]
    lower=find(frequencies<=target,1,'last'); upper=find(frequencies>target,1);
    selected=[selected lower upper];
end
selected=unique(selected);
end

function writeShapes(destination,snapshot)
selected=selectedModes(snapshot.frequencies);
fid=fopen(fullfile(destination,sprintf('shapes_n%03d.csv',snapshot.n)),'wt'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'mode,frequency_hz,member,x_m,y_m,ux_normalized,uy_normalized\n');
for k=selected
    u=reshape(snapshot.shapes(:,k),2,[]).'; u=u/max(sqrt(sum(u.^2,2)));
    for j=1:size(u,1)
        fprintf(fid,'%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g\n',k,snapshot.frequencies(k),1+(j>201),snapshot.points(j,:),u(j,:));
    end
end
end

function makePlots(destination,analytical,frames)
f=figure('visible','off'); cleanup=onCleanup(@() close(f));
for kind=0:1
    clf(f); a=analytical(analytical(:,1)==kind,:);
    semilogy(a(:,2),a(:,3:4),'o-'); grid on;
    xlabel('Elements'); ylabel('Absolute relative frequency error'); legend('Mode 1','Mode 2');
    names={'axial','bending'}; title([names{kind+1} ' analytical convergence']);
    print(f,fullfile(destination,[names{kind+1} '_convergence.svg']),'-dsvg');
end
clf(f); ns=cellfun(@(x)x.n,frames); first=cellfun(@(x)x.frequencies(1),frames); second=cellfun(@(x)x.frequencies(2),frames);
plot(ns,[first(:) second(:)],'o-'); grid on; xlabel('Elements per member'); ylabel('Frequency (Hz)'); legend('Mode 1','Mode 2');
print(f,fullfile(destination,'frame_convergence.svg'),'-dsvg');
last=frames{end}; clf(f); stem(last.frequencies,ones(size(last.frequencies)),'filled');
xlim([0 1700]); xlabel('Frequency (Hz)'); ylabel('Mode'); title(sprintf('Frame n=%d: %d modes below 1600 Hz',last.n,sum(last.frequencies<=1600))); grid on;
print(f,fullfile(destination,'spectrum.svg'),'-dsvg');
for k=selectedModes(last.frequencies)
    clf(f); hold on;
    u=reshape(last.shapes(:,k),2,[]).'; u=.12*u/max(sqrt(sum(u.^2,2)));
    for member=1:2
        ids=(member-1)*201+(1:201); p=last.points(ids,:);
        plot(p(:,1),p(:,2),'k--'); plot(p(:,1)+u(ids,1),p(:,2)+u(ids,2),'b-');
    end
    axis equal; grid on; xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Mode %d: %.6g Hz (amplitude scaled)',k,last.frequencies(k)));
    print(f,fullfile(destination,sprintf('frame_mode_%03d.svg',k)),'-dsvg');
end
end
