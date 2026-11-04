function vector = spatialNodalLoad(model,record)
%SPATIALNODALLOAD Conjugate projection; reject work in absent rotations.
node=record(2);
if node~=fix(node) || node<1 || node>model.numberOfNodes || numel(record)<8 || any(~isfinite(record))
    error('MKEF:InvalidNodalLoad','Invalid spatial nodal load.');
end
G=model.nodalKinematics((node-1)*6+(1:6),:); f=record(3:8).';
projection=G*(G.'*f);
if norm(projection(4:6)-f(4:6))>1e-10*norm(f(4:6))
    error('MKEF:InactiveRotationLoad','Node %d: moment acts in an inactive rotation direction.',node);
end
vector=G.'*f;
end
