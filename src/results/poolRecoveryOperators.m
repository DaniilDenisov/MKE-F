function model = poolRecoveryOperators(model, storage)
%POOLRECOVERYOPERATORS Share exactly equal local K/M pairs, never rounded values.
% Inline is retained when a table would enlarge compact JSON. IDs follow the
% first occurrence in element order. Hash matches are verified with isequal.
if strcmp(storage,'inline'),return;end
elements=model.elements; count=numel(elements); keys=cell(count,1);
if count<2 && strcmp(storage,'auto'),return;end
for i=1:count
    op=elements(i).recoveryOperators;
    bits=num2hex([op.K(:);op.M(:)]);
    keys{i}=hash('sha256',reshape(bits.',1,[]));
end
[~,first,group]=unique(keys);
% With no repeats a table only adds IDs and references. Avoid constructing
% and serializing it on the common all-unique path.
if numel(first)==count && strcmp(storage,'auto'),return;end
[first,order]=sort(first); inverse=zeros(numel(order),1);inverse(order)=1:numel(order);
ids=inverse(group); table=struct('id',{},'K',{},'M',{}); lengths=zeros(numel(first),1);
for j=1:numel(first)
    op=elements(first(j)).recoveryOperators;
    table(j)=struct('id',j,'K',op.K,'M',op.M);
    lengths(j)=numel(encodeResultJson(op));
end
for i=1:count
    op=elements(i).recoveryOperators;representative=table(ids(i));
    if ~isequal(op.K,representative.K) || ~isequal(op.M,representative.M)
        % Never merge unequal matrices, even in the event of a hash collision.
        error('MKEF:PostprocessorInvalidModel','Recovery operator hash collision; use inline storage.');
    end
end
if strcmp(storage,'auto')
    inlineBytes=sum(lengths(ids))+count*numel(',"recoveryOperators":');
    entries=arrayfun(@(op) op,table,'UniformOutput',false);
    tableBytes=numel(',"recoveryOperators":')+numel(encodeResultJson(entries));
    for i=1:count
        tableBytes=tableBytes+numel(',"recoveryOperatorId":')+numel(sprintf('%d',ids(i)));
    end
    if tableBytes>=inlineBytes,return;end
end
model.elements=rmfield(elements,'recoveryOperators');
for i=1:count,model.elements(i).recoveryOperatorId=ids(i);end
model.recovery.version=2;
model.recoveryOperators=table;
end
