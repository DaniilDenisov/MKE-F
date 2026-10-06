function manifest = generate_nafems_challenge5_cases(destination)
%GENERATE_NAFEMS_CHALLENGE5_CASES Deterministic, standard-format benchmark inputs.
% Explicit destination required; the study runner never calls this generator.
if nargin ~= 1 || ~ischar(destination) || isempty(destination)
    error('NAFEMS5:Destination', 'Supply an explicit output directory.');
end
if exist(destination, 'dir') ~= 7, mkdir(destination); end
manifest = struct('file', {}, 'kind', {}, 'n', {}, 'profile', {});
kinds = {'axial', 'bending', 'frame', 'frame_static'};
counts = {[1 2 4 8 16], [2 4 8 16], [1 2 4 8 16 25 50 100 200 400], 2};
for family = 1:numel(kinds)
    kind = kinds{family};
    for n = counts{family}
        name = sprintf('%s_n%03d.txt', kind, n);
        profile = 'quick';
        if n > 16, profile = 'full'; end
        if n == 400, profile = 'conditional'; end
        writeCase(fullfile(destination, name), kind, n);
        manifest(end+1) = struct('file', name, 'kind', kind, 'n', n, 'profile', profile);
    end
end
fid = fopen(fullfile(destination, 'manifest.csv'), 'wb');
assert(fid >= 0); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'file,kind,elements_per_member,profile\n');
for row = manifest
    fprintf(fid, '%s,%s,%d,%s\n', row.file, row.kind, row.n, row.profile);
end
end

function writeCase(path, kind, n)
E = 210e9; rho = 7800; radius = .01;
A = pi*radius^2; I = pi*radius^4/4;
isFrame = strncmp(kind, 'frame', 5);
isStatic = strcmp(kind, 'frame_static');
if isFrame
    nodes = [0 0; 0 1; 1 1]; elements = zeros(2*n, 2);
    ends = [2 3; 1 3];
    for member = 1:2
        chain = ends(member, 1);
        for k = 1:n-1
            nodes(end+1,:) = (1-k/n)*nodes(ends(member,1),:) + (k/n)*nodes(3,:);
            chain(end+1) = size(nodes,1);
        end
        chain(end+1) = 3;
        elements((member-1)*n+(1:n),:) = [chain(1:end-1).' chain(2:end).'];
    end
    supports = [4 1; 4 2];
    releases = [1 1; n 2; n+1 1; 2*n 2];
else
    nodes = [(0:n).'/n zeros(n+1,1)];
    elements = [(1:n).' (2:n+1).'];
    if strcmp(kind, 'axial')
        supports = [1 1; 2*ones(n,1) (2:n+1).'];
        releases = zeros(0,2);
    else
        supports = [4 1; 4 n+1]; releases = [1 1; n 2];
    end
end
fid = fopen(path, 'wb'); assert(fid >= 0); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# NAFEMS Challenge 5; SI: m, kg, s, N, Pa; frame element 113\n');
fprintf(fid, '# %s; elements per member: %d; no rotary rho*I mass term\n', kind, n);
analysis = 'modal'; if isStatic, analysis = 'static'; end
fprintf(fid, 'analysis\n%s\nnodes\n%d\n', analysis, size(nodes,1));
fprintf(fid, '%.17g,%.17g,0\n', nodes.');
fprintf(fid, 'elems_113\n%d\n', size(elements,1));
for k = 1:size(elements,1)
    fprintf(fid, '113,%d,%d,%.17g,%.17g,%.17g,%.17g\n', elements(k,:), A,E,rho,I);
end
if ~isempty(releases)
    fprintf(fid, 'releases\n%d\n', size(releases,1));
    fprintf(fid, '%d,%d,Mz\n', releases.');
end
fprintf(fid, 'bcfix\n%d\n', size(supports,1));
fprintf(fid, '%d,%d,0,0,0\n', supports.');
if isStatic, fprintf(fid, 'bcforce_stat\n1\n10,3,0,-1,0\n'); end
end
