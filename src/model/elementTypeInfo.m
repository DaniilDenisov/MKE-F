function info = elementTypeInfo(code)
%ELEMENTTYPEINFO Shared explicit registry. Never infer semantics from digits.
persistent catalog;
if isempty(catalog)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    catalog = jsondecode(fileread(fullfile(root,'shared','element-types.json')));
    if catalog.version ~= 1 || numel(unique([catalog.types.code])) ~= numel(catalog.types)
        error('MKEF:InvalidElementRegistry','Invalid element registry.');
    end
end
if nargin == 0, info = catalog.types; return; end
if ~isnumeric(code) || ~isscalar(code) || ~isreal(code) || ~isfinite(code)
    error('MKEF:UnsupportedElementType','Element code must be a finite scalar.');
end
index = find([catalog.types.code] == code);
if isempty(index), error('MKEF:UnsupportedElementType','Unsupported structural element type %g.',code); end
info = catalog.types(index);
info.dofPerNode = numel(info.dofLabels);
end
