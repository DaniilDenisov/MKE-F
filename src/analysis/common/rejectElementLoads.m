function rejectElementLoads(model)
%REJECTELEMENTLOADS Prevent silently ignoring static-only element loads.
if isfield(model, 'elementLoads') && ~isempty(model.elementLoads)
    error('MKEF:AnalysisLoadMismatch', 'Element loads are supported only in static analysis.');
end
end
