% Run format: octave-cli */run_svd_task.m input_path input_name output_path
arg_list=argv();
mat_name=[arg_list{1},arg_list{2}];
A=open(mat_name);
A=A.data;
label=label.data;
[batch,m,n]=size(A);
u=zeros(batch,m,m);
v=zeros(batch,n,n);
s=zeros(batch,m,n);
run_time=linspace(0,0,200);
size(run_time)
for i=1:200
    a=squeeze(A(i,:,:));
    tic
    [u_t,s_t,v_t]=my_svd_2(a);
    run_time(i)=toc;
    u(i,:,:)=u_t;
    s(i,:,:)=s_t;
    v(i,:,:)=v_t;
end
disp(run_time);
res_name=[arg_list{3},'res.mat'];
save(res_name,'u','v','s','run_time');
