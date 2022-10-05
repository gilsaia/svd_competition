% Run format: octave-cli */run_svd_task.m input_path input_name label_name output_path
arg_list=argv();
mat_name=[arg_list{1},arg_list{2}];
A=open(mat_name);
A=A.data;
label_name=[arg_list{1},arg_list{3}];
label=open(label_name);
label=label.data;
[batch,m,n]=size(A);
E1=zeros(1,batch);
G1=zeros(1,batch);
run_time=linspace(0,0,200);
for i=1:200
    a=squeeze(A(i,:,:));
    l=squeeze(label(i,:));
    tic
    [u_t,s_t,v_t]=my_svd_1(a,l);
    run_time(i)=toc;
    [e1,g1]=check_svd(a,u_t,s_t,v_t);
    E1(i)=e1;
    G1(i)=g1;
    disp(sprintf('Index:%d Rank:%d E1:%f G1:%f Time:%f',i,l,e1,g1,run_time(i)));
end
res_name=[arg_list{4},'res.mat'];
save('-v6',res_name,'E1','G1','run_time');
