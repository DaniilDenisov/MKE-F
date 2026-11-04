function validateSpatialFrameSections(model)
%VALIDATESPATIALFRAMESECTIONS Assignment policy for imported models and all calculation paths.
if ~isfield(model,'elementData') || ~isfield(model,'sections'), return; end
allowed={'rectangle','circle','square_tube','rectangular_tube','circular_tube'};
members={};
for i=1:numel(model.elementData)
    e=model.elementData(i);
    if e.type~=1436 || ~isfield(e,'sectionId') || e.sectionId==0, continue; end
    index=find([model.sections.id]==e.sectionId,1);
    if isempty(index), continue; end % Catalog consistency is checked when building the model.
    s=model.sections(index);
    if ~ismember(s.type,allowed), members{end+1}=sprintf('%d (%s)',i,s.type); end
end
if ~isempty(members)
    error('MKEF:UnsupportedFrameSection', ...
        ['Unsupported sections on members %s. Frames 1436 support only rectangles, ' ...
         'solid circles and square, rectangular or circular tubes. Open profiles require ' ...
         'shear-center offsets and restrained warping that this model does not include.'], ...
        strjoin(members,', '));
end
end
