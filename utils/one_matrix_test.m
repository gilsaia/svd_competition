% Run format: octave-cli */one_matrix_test.m task
arg_list=argv();
[A]=generate_matrix();
if arg_list{1}=='1'
    [U,S,V]=my_svd_1(A,32);
elseif arg_list{1}=='2'
    [U,S,V]=my_svd_2(A);
elseif arg_list{1}=='3'
    [U,S,V]=my_svd_3(A);
elseif arg_list{1}=='4'
    [inv]=my_inverse(A);
end
disp('Run Complete');