import argparse
import os
import numpy as np
from scipy.io import loadmat


sh_dict = {1: 'utils/run_svd_task_with_r.m',
           2: 'utils/run_svd_task.m',
           3: 'utils/run_svd_task_dense.m',
           4: 'utils/run_inv_task.m'}

e1_dict = {1: 1e-7, 2: 1e-7, 3: 1e-2}
g1_dict = {1: 1e-4, 2: 1e-4, 3: 1e-4}
e2_dict = {4: 1e-6}


def check_svd(args, truth, u, s, v):
    e1_thd = e1_dict[args.task]
    g1_thd = g1_dict[args.task]

    vh = v.T.conj()
    res = np.dot(u, np.dot(s, vh))
    remain = truth-res
    remain_norm = np.linalg.norm(remain, 'fro')
    truth_norm = np.linalg.norm(truth, 'fro')
    e1 = remain_norm/truth_norm

    uh = u.T.conj()
    ut = np.dot(uh, u)
    iu = np.eye(ut.shape[0])
    g1u = np.linalg.norm(iu-ut, 'fro')/np.linalg.norm(iu, 'fro')
    vt = np.dot(vh, v)
    iv = np.eye(vt.shape[0])
    g1v = np.linalg.norm(iv-vt, 'fro')/np.linalg.norm(iv, 'fro')
    g1 = max(g1u, g1v)

    if e1 > e1_thd or g1 > g1_thd:
        return False, e1, g1
    return True, e1, g1


def check_inv(args):
    raise NotImplementedError()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str, required=True)
    parser.add_argument('--output_path', type=str, required=True)
    parser.add_argument('--save_path', type=str, default='runs/')
    parser.add_argument('--project_name', type=str, default='svd_competition')
    parser.add_argument('--entity_name', type=str, default='svd_competition')
    parser.add_argument('--run_name', type=str, default='origin')
    parser.add_argument('--simple_check', action='store_true')
    parser.add_argument('--complete_check', action='store_true')
    parser.add_argument('--measure', action='store_true')
    parser.add_argument('--task', type=int, choices=[1, 2, 3, 4], default=1)
    return parser.parse_args()


def prepare_dir(args):
    index = 0
    run_name = f'{args.run_name}-{index}'
    run_path = f'{args.save_path}{run_name}'
    while os.path.exists(run_path):
        index += 1
        run_name = f'{args.run_name}-{index}'
        run_path = f'{args.save_path}{run_name}'
    if args.measure:
        os.makedirs(run_path)
    args.run_name = run_name

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)


def measure(args):
    raise NotImplementedError()


def complete_check(args):
    raise NotImplementedError()


def simple_check(args):
    data_name = 'data_m_256_n_128_dataNum_200.mat'
    label_name = 'data_m_256_n_128_label_r.mat'
    sh = sh_dict[args.task]
    input_name = f'{args.input_path}{data_name}'
    output_name = f'{args.output_path}res.mat'
    if args.task == 1:
        cmd_args = f'{args.input_path} {data_name} {label_name} {args.output_path}'
    else:
        cmd_args = f'{args.input_path} {data_name} {args.output_path}'
    cmd = f'octave-cli --path code/ {sh} {cmd_args}'
    os.system(cmd)
    print('Run algorithm end!')
    if not os.path.exists(output_name):
        print('Not find target file\nCheck error')
        return
    truth_mat = loadmat(input_name)
    truth_mat = truth_mat['data']
    output_res = loadmat(output_name)
    u = output_res['u']
    v = output_res['v']
    s = output_res['s']
    e1_sum = 0
    g1_sum = 0
    e2_sum = 0
    for i in range(truth_mat.shape[0]):
        if args.task == 4:
            check, e2 = check_inv(args)
        else:
            check, e1, g1 = check_svd(args, truth_mat[i], u[i], s[i], v[i])
            if not check:
                print(
                    f'Find result exceed threshold\nIndex:{i}\tE1:{e1}\tG1:{g1}\nCheck error')
                return
            e1_sum += e1
            g1_sum += g1
    e1_avg = e1_sum/truth_mat.shape[0]
    g1_avg = g1_sum/truth_mat.shape[0]
    e2_avg = e2_sum/truth_mat.shape[0]
    print(
        f'Simple check pass!\nTask:{args.task}\tE1:{e1_avg}\tG1:{g1_avg}\tE2:{e2_avg}')


if __name__ == '__main__':
    args = get_args()
    prepare_dir(args)
    if args.measure:
        measure(args)
    elif args.complete_check:
        complete_check(args)
    elif args.simple_check:
        simple_check(args)
