function test_postprocessor_site()
%TEST_POSTPROCESSOR_SITE Verify offline assets and security invariants.

testsDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(testsDir);
postDir = fullfile(rootDir, 'postprocessor');
pages = {fullfile(postDir, 'index.html'), ...
    fullfile(postDir, 'tests', 'index.html')};

for pageIndex = 1:numel(pages)
    page = pages{pageIndex};
    assert(exist(page, 'file') == 2, ...
        ['MKEF:VerificationFailed: missing postprocessor page ' page '.']);
    html = fileread(page);
    assert(isempty(regexpi(html, '<script[^>]+type=["'']module', 'once')), ...
        'MKEF:VerificationFailed: postprocessor uses file-restricted modules.');
    assert(isempty(regexpi(html, ...
        '(?:src|href)\s*=\s*["''](?:https?:)?//', 'once')), ...
        'MKEF:VerificationFailed: postprocessor has an external runtime asset.');
    targets = regexp(html, '(?:src|href)\s*=\s*["'']([^"'']+)["'']', 'tokens');
    for i = 1:numel(targets)
        target = targets{i}{1};
        if isempty(target) || target(1) == '#'
            continue;
        end
        targetPath = fullfile(fileparts(page), strrep(target, '/', filesep));
        assert(exist(targetPath, 'file') == 2, ...
            ['MKEF:VerificationFailed: unresolved postprocessor asset ' target '.']);
    end
end

mainHtml = fileread(fullfile(postDir, 'index.html'));
preprocessorRenderer = fileread(fullfile(rootDir, 'preprocessor', 'js', 'spatial-panel.js'));
postprocessorRenderer = fileread(fullfile(postDir, 'js', 'webgl-renderer.js'));
assert(~isempty(strfind(mainHtml, '../shared/vendor/spatial.js')) && ...
       ~isempty(strfind(preprocessorRenderer, 'MKEFSpatial.SpatialScene')) && ...
       ~isempty(strfind(postprocessorRenderer, 'MKEFSpatial.SpatialScene')), ...
    'MKEF:VerificationFailed: processors do not share the WebGL scene.');
modalControlIds = {'mode-number', 'previous-mode', 'modal-play-pause', ...
    'next-mode', 'modal-phase', 'modal-phase-value', ...
    'modal-playback-speed', 'modal-amplitude', 'mode-frequency'};
for i = 1:numel(modalControlIds)
    assert(~isempty(strfind(mainHtml, ['id="' modalControlIds{i} '"'])), ...
        ['MKEF:VerificationFailed: missing modal control ' ...
         modalControlIds{i} '.']);
end

transientControlIds = {'time-index', 'previous-frame', 'play-pause', ...
    'next-frame', 'playback-speed', 'frame-stride', 'time-value', ...
    'transient-metadata', 'history-node', 'history-dof', ...
    'history-quantity', 'history-value'};
for i = 1:numel(transientControlIds)
    assert(~isempty(strfind(mainHtml, ['id="' transientControlIds{i} '"'])), ...
        ['MKEF:VerificationFailed: missing transient control ' ...
         transientControlIds{i} '.']);
end

assets = [dir(fullfile(postDir, 'js', '*.js')); ...
          dir(fullfile(postDir, 'assets', '*.css'))];
for i = 1:numel(assets)
    content = fileread(fullfile(assets(i).folder, assets(i).name));
    networkCheck = strrep(content, 'http://www.w3.org/2000/svg', '');
    if strcmp(assets(i).name,'member-maxima.js')
        % Optional explicit calculation uses the same-origin Compose service.
        % Local result viewing remains available without it (browser-tested).
        assert(~isempty(strfind(content,'const api=''/api/v1/sections/''')));
        assert(~isempty(strfind(content,'location.protocol===''file:''')));
        networkCheck=strrep(networkCheck,'fetch(api+','sectionRequest(');
    end
    assert(isempty(regexpi(networkCheck, 'https?://|@import\s|fetch\s*\(', 'once')), ...
        ['MKEF:VerificationFailed: external network dependency in ' assets(i).name '.']);
    if strcmp(assets(i).name(end-2:end), '.js')
        assert(isempty(regexp(content, ...
            'innerHTML|\beval\s*\(|new\s+Function\s*\(', 'once')), ...
            ['MKEF:VerificationFailed: unsafe DOM/code API in ' assets(i).name '.']);
    end
end

assert(exist(fullfile(postDir, 'schema-v1.md'), 'file') == 2, ...
    'MKEF:VerificationFailed: postprocessor schema documentation is missing.');
end
