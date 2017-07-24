function st=mergeKITTI(st1,st2)



[F1,N1]=size(st1.Xi);
[F2,N2]=size(st2.Xi);

assert(F1==F2,'states must be equal length');
assert(all(st1.frameNums==st2.frameNums),'frameNums must be identical');


% assume st1=car, st2=ped
st1.C=ones(F1,N1);st2.C=2*ones(F2,N2);

st=st1;

if isfield(st1,'Xi') && isfield(st2,'Xi')
    st.Xi=[st1.Xi st2.Xi];
    st.Yi=[st1.Yi st2.Yi];
    st.W=[st1.W st2.W];
    st.H=[st1.H st2.H];
end
if isfield(st1,'Y') && isfield(st2,'X')
    st.X=[st1.X st2.X];
    st.Y=[st1.Y st2.Y];
end

st.C=[st1.C st2.C];

st.sc1=st1.sceneInfo;st.sc2=st2.sceneInfo;
st.opt1=st1.opt;st.opt2=st2.opt;


%% merge different classes into one

% file1='/home/amilan/research/projects/dctracking/results/1202Kc-5/0000.txt';
% file2='/home/amilan/research/projects/dctracking/results/1204Kb-2/0000.txt';



