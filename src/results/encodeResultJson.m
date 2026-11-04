function json = encodeResultJson(value, pretty)
%ENCODERESULTJSON Finite IEEE doubles with round-trip precision, explicit arrays.
% Octave 11.3 jsonencode rounds sub-epsilon doubles to zero. Use it only for
% strings/logicals; %.17g preserves nonzero doubles including subnormals.
% Normalize signed zero, as JSON.stringify does on browser download.
% The caller prepares singleton schema arrays as cells before encoding.
if nargin<2, pretty=false; end
json=encode(value,0,pretty);
end

function json=encode(v,level,pretty)
if isnumeric(v)
    if ~isreal(v) || issparse(v) || any(~isfinite(v(:))) || ndims(v)>2
        error('MKEF:PostprocessorEncodingFailed','Expected finite real dense numbers.');
    end
    if isempty(v),json='[]';
    elseif isscalar(v),if v==0,json='0';else,json=sprintf('%.17g',v);end
    elseif isvector(v),json=numbers(v);
    else
        rows=cell(size(v,1),1);
        for i=1:size(v,1),rows{i}=numbers(v(i,:));end
        json=['[' strjoin(rows,',') ']'];
    end
elseif islogical(v) || ischar(v)
    json=jsonencode(v);
elseif iscell(v)
    if isempty(v),json='[]';return;end
    % Prepared numerical rows are cells. Format them together, not recursively
    % one double at a time; this keeps long transient histories practical.
    if all(cellfun('isclass',v(:),'double') & cellfun('numel',v(:))==1)
        a=[v{:}];
        if ~isreal(a)||any(~isfinite(a)),error('MKEF:PostprocessorEncodingFailed','Nonfinite numeric row.');end
        json=numbers(a);
    else
        parts=cellfun(@(x) encode(x,level+1,pretty),v(:),'UniformOutput',false);
        json=container(parts,'[',']',level,pretty);
    end
elseif isstruct(v)
    if numel(v)~=1
        json=encode(arrayfun(@(x) x,v,'UniformOutput',false),level,pretty);return;
    end
    fields=fieldnames(v);parts=cell(numel(fields),1);separator=':';if pretty,separator=': ';end
    for i=1:numel(fields)
        parts{i}=[jsonencode(fields{i}) separator encode(v.(fields{i}),level+1,pretty)];
    end
    json=container(parts,'{','}',level,pretty);
else
    error('MKEF:PostprocessorEncodingFailed','Unsupported JSON value type %s.',class(v));
end
end
function json=numbers(v)
v(v==0)=0;
text=sprintf('%.17g,',v);json=['[' text(1:end-1) ']'];
end
function json=container(parts,left,right,level,pretty)
if pretty && ~isempty(parts)
    pad=repmat(' ',1,2*(level+1));sep=sprintf(',\n');
    json=[left sprintf('\n') pad strjoin(parts,[sep pad]) sprintf('\n') repmat(' ',1,2*level) right];
else,json=[left strjoin(parts,',') right];end
end
