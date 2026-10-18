function test_sections()
% Independent profile values and beam solutions; fixtures shared with JS.
types={'rectangle','circle','channel','tee','ibeam','square_tube','rectangular_tube','circular_tube'};
output=struct('section',{},'angle',{},'properties',{},'stress',{});
for k=1:numel(types)
    s=struct('id',k,'name',types{k},'type',types{k},'dimensions',struct('b',.1,'h',.2,'tw',.008,'tf',.012,'d',.1,'t',.005));
    for angle=[0 90 180 270]
        p=sectionProperties(s,angle);
        if strcmp(s.type,'rectangle')
            near(p.area,.02);near(p.Iy,.2*.1^3/12);near(p.Iz,.1*.2^3/12);
        elseif strcmp(s.type,'circle')
            near(p.area,pi*.1^2/4);near(p.momentOfInertia,pi*.1^4/64);
        elseif strcmp(s.type,'circular_tube')
            near(p.area,pi*(.1^2-.09^2)/4);near(p.momentOfInertia,pi*(.1^4-.09^4)/64);
            near(p.centroidY,0);near(p.centroidZ,0);near(p.etaMin,-.05);near(p.etaMax,.05);
        elseif ismember(s.type,{'square_tube','rectangular_tube'})
            h=.2; if strcmp(s.type,'square_tube'), h=.1; end
            near(p.area,.1*h-.09*(h-.01));
            near(p.Iy,(h*.1^3-(h-.01)*.09^3)/12);
            near(p.Iz,(.1*h^3-.09*(h-.01)^3)/12);
            near(p.centroidY,h/2);near(p.centroidZ,.05);
            extent=h/2; expectedI=p.Iz;
            if ismember(angle,[90 270]), extent=.05; expectedI=p.Iy; end
            near(p.etaMin,-extent);near(p.etaMax,extent);near(p.momentOfInertia,expectedI);
        elseif strcmp(s.type,'tee')
            near(p.area,.1*.012+.008*.188);
            near(p.centroidY,(.1*.012*.194+.008*.188*.094)/p.area);
            near(p.Iy,(.012*.1^3+.188*.008^3)/12);
        else
            near(p.area,2*.1*.012+.008*.176);
            near(p.Iz,(.1*.2^3-(.1-.008)*.176^3)/12);
            expectedIy=(2*.012*.1^3+.176*.008^3)/12;
            if strcmp(s.type,'channel')
                cz=(2*.1*.012*.05+.008*.176*.004)/p.area;
                near(p.centroidZ,cz);
                expectedIy=expectedIy+2*.1*.012*(.05-cz)^2+.008*.176*(.004-cz)^2;
            end
            near(p.Iy,expectedIy);
        end
        stress=recoverSectionStress(p,113,2,[-40;150;-30;0;0;0],[20 -300 100 200]);
        output(end+1)=struct('section',s,'angle',angle,'properties',p,'stress',stress);
        % A dense scan independently bounds the analytic global extrema.
        for xi=linspace(0,1,301)
            x=2*xi; N=40-20*x-40*x^2/2; M=30+150*x-300*x^2/2+250*x^3/6;
            values=N/p.area-M*[p.etaMin p.etaMax]/p.momentOfInertia;
            assert(min(values)>=stress.minimum.value-1e-6 && max(values)<=stress.maximum.value+1e-6);
        end
    end
end
if ~exist('output','dir'), mkdir('output'); end
% Near-endpoint moment roots must not replace exact normalized endpoints.
profile=sectionProperties(struct('type','rectangle','dimensions',struct('b',.01,'h',.02)));
for root=[1e-12 1-1e-12 0 1]
    endpointStress=recoverSectionStress(profile,113,1,[0;1;root;0;0;0],[0 0 0 0]);
    assert(isequal(endpointStress.positions,[0 1]));
end
fixture=fullfile(fileparts(mfilename('fullpath')),'fixtures','CaseSectionEndpointRounding.txt');
endpointProblem=StructFEProblem(fixture,struct('verbose',false,'plotting',false));
endpointResult=endpointProblem.RunStatic();
for k=1:numel(endpointResult.elementResults)
    positions=endpointResult.elementResults(k).sectionStress.positions;
    assert(positions(1)==0 && positions(end)==1 && all(diff(positions)>0));
end
exportPostprocessorData(endpointProblem.GetAnalysisModel(),endpointResult,fullfile('output','SectionCaseEndpointRounding.json'));
fid=fopen(fullfile('output','section-parity.json'),'w');fprintf(fid,'%s',jsonencode(output));fclose(fid);
p=StructFEProblem(exampleCasePath('CaseSectionBeam.txt'),struct('verbose',false,'plotting',false));
m=p.GetAnalysisModel();r=p.RunStatic(); st=r.elementResults.sectionStress;
near(st.maximum.value,750000); near(st.minimum.value,-750000); near(st.maximum.xi,.5);
exportPostprocessorData(m,r,fullfile('output','CaseSectionBeam.json'));
text=fileread(exampleCasePath('CaseSectionBeam.txt'));
% Tube geometry drives beam stiffness, mass, stress and exported results.
for type={'square_tube','rectangular_tube','circular_tube'}
    section=struct('id',1,'name',type{1},'type',type{1},'dimensions',struct('b',.1,'h',.2,'d',.1,'t',.005));
    source=regexprep(text,'\{"id":1[^\r\n]+',jsonencode(section));
    for angle=[0 90 180 270]
        tube=fromText(strrep(source,sprintf('1,1,0\n'),sprintf('1,1,%d\n',angle)));
        tubeResult=tube.RunStatic(); tubeModel=tube.GetAnalysisModel();
        props=sectionProperties(section,angle);
        near(tubeResult.displacements(3),r.displacements(3)*(.1*.2^3/12)/props.momentOfInertia);
        near(full(tube.M),full(p.M)*props.area/.02);
        near(tubeResult.elementResults.sectionStress.maximum.value,500*props.etaMax/props.momentOfInertia);
        exportPostprocessorData(tubeModel,tubeResult,fullfile('output',['SectionCase_' type{1} '_' num2str(angle) '.json']));
    end
    for thickness=[0 -.01 .05 .1 NaN Inf]
        bad=section;bad.dimensions.t=thickness;mustFail(@() sectionProperties(bad));
    end
    bad=section;bad.dimensions=rmfield(bad.dimensions,'t');mustFail(@() sectionProperties(bad));
end
bad=struct('type','rectangular_tube','dimensions',struct('b',.2,'h',.1,'t',.05));
mustFail(@() sectionProperties(bad));
% Manual A/I produce the same solution with an unused catalog or no catalog.
manual=regexprep(text,'element_sections\s+1\s+1,1,0\s*','');
for withCatalog=[true false]
    source=manual;
    if ~withCatalog, source=regexprep(source,'sections\s+1\s+[^\r\n]+\s*',''); end
    direct=fromText(source); directResult=direct.RunStatic(); directModel=direct.GetAnalysisModel();
    assert(directModel.elementData(1).sectionId==0);
    near(full(direct.K),full(p.K)); near(full(direct.M),full(p.M));
    near(directResult.displacements,r.displacements);
    near(directResult.elementResults.localEndForces,r.elementResults.localEndForces);
end
% Rotate the rectangular beam: weak-axis deflection is four times larger.
p90=fromText(strrep(text,sprintf('1,1,0\n'),sprintf('1,1,90\n')));r90=p90.RunStatic();
near(r90.displacements(3),4*r.displacements(3));
near(r90.elementResults.sectionStress.maximum.value,2*st.maximum.value);
% Explicit properties are overridden by geometry, including mass.
pOther=fromText(strrep(text,'0.02,200000000000,7800,0.00006666666666666667','1,200000000000,7800,1'));
near(full(pOther.K),full(p.K));near(full(pOther.M),full(p.M));
% V6 must compose with every existing schema feature and analysis mode.
for name={'CaseLinearFrame','CaseLinearRelease','CaseLinearMPC','CaseReleaseModal','CaseReleaseTransient','CaseMPCModal','CaseMPCTransient'}
    source=fileread(exampleCasePath([name{1} '.txt']));
    source=[source sprintf('\nsections\n1\n{"id":1,"name":"T","type":"tee","dimensions":{"b":0.1,"h":0.2,"tw":0.008,"tf":0.012}}\nelement_sections\n1\n1,1,270\n')];
    t=fromText(source);model=t.GetAnalysisModel();
    result=t.RunSelected();
    data=exportPostprocessorData(model,result,fullfile('output',['Section' name{1} '.json']));assert(data.version==6);
end
bad=struct('type','tee','dimensions',struct('b',.1,'h',.2,'tw',.1,'tf',.01));mustFail(@() sectionProperties(bad));
mustFail(@() fromText(strrep(text,'1,1,0','1,1,45')));
mustFail(@() fromText([text sprintf('\nelement_sections\n1\n1,1,0\n')]));
mustFail(@() fromText(strrep(text,'1,1,0','1,99,0')));
end
function p=fromText(text)
file=[tempname() '.txt'];fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);cleanup=onCleanup(@() delete(file));
p=StructFEProblem(file,struct('verbose',false,'plotting',false));
end
function near(a,b)
assert(norm(a(:)-b(:),inf)<=1e-9*max(norm(b(:),inf),1e-15));
end
function mustFail(action)
try,action();catch,return;end
error('Expected failure.');
end
