% Run format: octave-cli */run_inv_task.m input_path input_name output_path
arg_list=argv();
mat_name=[arg_list{1},arg_list{2}];
A=open(mat_name);
A=A.data;
[batch,m,n]=size(A);
inv=zeros(batch,m,m);
run_time=linspace(0,0,200);
for i=1:200
    a=squeeze(A(i,:,:));
    tic
    [inv_a]=my_inverse(a);
    run_time(i)=toc;
    inv(i,:,:)=inv_a;
end
res_name=[arg_list{3},'res.mat'];
save('-v6',res_name,'inv','run_time');
