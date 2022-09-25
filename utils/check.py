import argparse
import os
import shutil
import numpy as np
import wandb
from scipy.io import loadmat

from utils import AverageMeter


sh_dict = {1: 'utils/run_svd_task_with_r.m',
           2: 'utils/run_svd_task.m',
           3: 'utils/run_svd_task_dense.m',
           4: 'utils/run_inv_task.m'}

e1_dict = {1: 1e-7, 2: 1e-7, 3: 1e-2}
g1_dict = {1: 1e-4, 2: 1e-4, 3: 1e-4}
e2_dict = {4: 1e-6}


def check_svd(args, truth, u, s, v):
    # deprecated
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


def check_inv(args, truth, inv):
    # deprecated
    e2_thd = e2_dict[args.task]

    temp = np.eye(truth.shape[0], truth.shape[0])
    temp = temp+np.dot(truth, truth.T.conj())
    res = np.dot(temp, inv)-np.eye(truth.shape[0], truth.shape[0])
    e2 = np.linalg.norm(res, 'fro')/np.sqrt(truth.shape[0])

    if e2 > e2_thd:
        return False, e2
    return True, e2


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
    parser.add_argument('--matlab', action='store_true')
    parser.add_argument('--task', type=int, choices=[1, 2, 3, 4], default=1)
    parser.add_argument('--wandb_mode', type=str, default='online')
    parser.add_argument('--data_len', type=int, default=3, choices=[1, 2, 3])
    return parser.parse_args()


DATA_LEN = 0  # don't change here, change args.data_len


def prepare(args):
    index = 0
    run_name = f'{args.run_name}-{index}'
    run_path = f'{args.save_path}{run_name}'
    code_sync_path = f'{args.save_path}{args.run_name}'
    while os.path.exists(run_path):
        index += 1
        run_name = f'{args.run_name}-{index}'
        run_path = f'{args.save_path}{run_name}'
    if args.measure:
        shutil.copytree('code/', code_sync_path, dirs_exist_ok=True)
        os.makedirs(run_path)
    args.run_name = run_name
    global DATA_LEN
    DATA_LEN = args.data_len

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)


def clean_tmpfile(args):
    os.remove(f'{args.output_path}res.mat')


data_name_list = ['data_m_256_n_128_dataNum_200.mat',
                  'data_m_512_n_256_dataNum_200.mat',
                  'data_m_1024_n_512_dataNum_200.mat']

label_name_list = ['data_m_256_n_128_label_r.mat',
                   'data_m_512_n_256_label_r.mat',
                   'data_m_1024_n_512_label_r.mat']

record_name_list = ['m_256_n_128', 'm_512_n_256', 'm_1024_n_512']


def get_run_cmd(args, data_name, label_name):
    sh = sh_dict[args.task]
    if args.task == 1:
        cmd_args = f'{args.input_path} {data_name} {label_name} {args.output_path}'
    else:
        cmd_args = f'{args.input_path} {data_name} {args.output_path}'
    if args.matlab:
        sh = sh.rstrip('.m')
        cmd = f'matlab -nodesktop -nosplash --path code/ --path utils/ -r {sh} {cmd_args}'
    else:
        cmd = f'octave-cli --path code/ --path utils/ {sh} {cmd_args}'
    return cmd


RETRY_NUM = 3
MATRIX_NUM = 200


def measure(args):
    save_path = f'{args.save_path}{args.run_name}'
    args.wandb_name = f'{args.run_name}-task{args.task}'

    if args.task != 4:
        summary_record_dict = {'e1': {}, 'g1': {}, 'time': {}}
    else:
        summary_record_dict = {'e2': {}, 'time': {}}
    for (key, item) in summary_record_dict.items():
        for j in range(DATA_LEN):
            record = record_name_list[j]
            item[record] = []

    for i in range(RETRY_NUM):
        job_name = f'{args.wandb_name}-measure-{i}'
        wandb.init(project=args.project_name,
                   entity=args.entity_name, config=args, group=args.wandb_name, job_type='measure', name=job_name, dir=save_path, mode=args.wandb_mode)
        if args.task != 4:
            wandb_record_dict = {'e1': {}, 'g1': {}, 'time': {}}
        else:
            wandb_record_dict = {'e2': {}, 'time': {}}
        for j in range(DATA_LEN):
            data_name = data_name_list[j]
            label_name = label_name_list[j]
            input_name = f'{args.input_path}{data_name}'
            output_name = f'{args.output_path}res.mat'
            if not os.path.exists(input_name):
                continue

            record_name = record_name_list[j]
            for (key, item) in wandb_record_dict.items():
                item[record_name] = []

            cmd = get_run_cmd(args, data_name, label_name)
            print(f'Run Command:{cmd}')
            os.system(cmd)
            print('Run Command end!')
            if not os.path.exists(output_name):
                print('Not find target file\nCheck error')
                return

            e1_meter = AverageMeter()
            e2_meter = AverageMeter()
            g1_meter = AverageMeter()
            time_meter = AverageMeter()
            if args.task != 4:
                truth_mat = loadmat(input_name)
                truth_mat = truth_mat['data']
                output_res = loadmat(output_name)
                E1 = np.squeeze(output_res['E1'])
                G1 = np.squeeze(output_res['G1'])
                runtime = np.squeeze(output_res['run_time'])
                e1_thd = e1_dict[args.task]
                g1_thd = g1_dict[args.task]
                for t in range(truth_mat.shape[0]):
                    e1 = E1[t]
                    g1 = G1[t]
                    check = e1 < e1_thd and g1 < g1_thd
                    if not check:
                        print(
                            f'Find result exceed threshold\nTask:{args.task}\tIndex:{t}\tE1:{e1}\tG1:{g1}\nCheck error')
                        return
                    e1_meter.update(e1)
                    g1_meter.update(g1)
                    time_meter.update(runtime[t])

                    wandb_record_dict['e1'][record_name].append(e1)
                    wandb_record_dict['g1'][record_name].append(g1)
                    wandb_record_dict['time'][record_name].append(runtime[t])

                summary_record_dict['e1'][record_name].append(e1_meter)
                summary_record_dict['g1'][record_name].append(g1_meter)
                summary_record_dict['time'][record_name].append(time_meter)
            else:
                truth_mat = loadmat(input_name)
                truth_mat = truth_mat['data']
                output_res = loadmat(output_name)
                E2 = np.squeeze(output_res['E2'])
                runtime = np.squeeze(output_res['run_time'])
                e2_thd = e2_dict[args.task]
                for t in range(truth_mat.shape[0]):
                    e2 = E2[t]
                    check = e2 < e2_thd
                    if not check:
                        print(
                            f'Find result exceed threshold\nTask:{args.task}\tIndex:{t}\tE2:{e2}\nCheck error')
                        return
                    e2_meter.update(e2)
                    time_meter.update(runtime[t])

                    wandb_record_dict['e2'][record_name].append(e2)
                    wandb_record_dict['time'][record_name].append(runtime[t])

                summary_record_dict['e2'][record_name].append(e2_meter)
                summary_record_dict['time'][record_name].append(time_meter)
            print(
                f'Complete check pass!\nTask:{args.task}\tData:{data_name}\tE1:{e1_meter.avg}\tG1:{g1_meter.avg}\tE2:{e2_meter.avg}\tRun time:{time_meter.avg}')

        for j in range(MATRIX_NUM):
            for (key, item) in wandb_record_dict.items():
                for (subkey, logs) in item.items():
                    wandb.log({f'{key}-{subkey}': logs[j]}, step=j)
        wandb.finish()

    wandb.init(project=args.project_name,
               entity=args.entity_name, config=args, group=args.wandb_name, job_type='summary', name=args.wandb_name, dir=save_path, mode=args.wandb_mode)

    for i in range(RETRY_NUM):
        for (key, item) in summary_record_dict.items():
            for (subkey, record) in item.items():
                wandb.log({f'Summary-{key}-{subkey}': record[i].avg}, step=i)

    for (key, item) in summary_record_dict.items():
        for (subkey, record) in item.items():
            sum = 0
            for i in range(RETRY_NUM):
                sum += record[i].avg
            sum /= RETRY_NUM
            wandb.run.summary[f'{key}-{subkey}'] = sum
    wandb.finish()


def complete_check(args):
    for i in range(DATA_LEN):
        data_name = data_name_list[i]
        label_name = label_name_list[i]
        input_name = f'{args.input_path}{data_name}'
        output_name = f'{args.output_path}res.mat'
        if not os.path.exists(input_name):
            continue
        cmd = get_run_cmd(args, data_name, label_name)
        print(f'Run Command:{cmd}')
        os.system(cmd)
        print('Run Command end!')
        if not os.path.exists(output_name):
            print('Not find target file\nCheck error')
            return
        e1_meter = AverageMeter()
        e2_meter = AverageMeter()
        g1_meter = AverageMeter()
        time_meter = AverageMeter()
        if args.task != 4:
            truth_mat = loadmat(input_name)
            truth_mat = truth_mat['data']
            output_res = loadmat(output_name)
            E1 = np.squeeze(output_res['E1'])
            G1 = np.squeeze(output_res['G1'])
            runtime = np.squeeze(output_res['run_time'])
            e1_thd = e1_dict[args.task]
            g1_thd = g1_dict[args.task]
            for j in range(truth_mat.shape[0]):
                e1 = E1[j]
                g1 = G1[j]
                check = e1 < e1_thd and g1 < g1_thd
                if not check:
                    print(
                        f'Find result exceed threshold\nTask:{args.task}\tIndex:{j}\tE1:{e1}\tG1:{g1}\nCheck error')
                    return
                e1_meter.update(e1)
                g1_meter.update(g1)
                time_meter.update(runtime[j])
        else:
            truth_mat = loadmat(input_name)
            truth_mat = truth_mat['data']
            output_res = loadmat(output_name)
            E2 = np.squeeze(output_res['E2'])
            runtime = np.squeeze(output_res['run_time'])
            e2_thd = e2_dict[args.task]
            for j in range(truth_mat.shape[0]):
                e2 = E2[j]
                check = e2 < e2_thd
                if not check:
                    print(
                        f'Find result exceed threshold\nTask:{args.task}\tIndex:{j}\tE2:{e2}\nCheck error')
                    return
                e2_meter.update(e2)
                time_meter.update(runtime[j])
        print(
            f'Complete check pass!\nTask:{args.task}\tData:{data_name}\tE1:{e1_meter.avg}\tG1:{g1_meter.avg}\tE2:{e2_meter.avg}\tRun time:{time_meter.avg}')


def simple_check(args):
    data_name = 'data_m_256_n_128_dataNum_200.mat'
    label_name = 'data_m_256_n_128_label_r.mat'
    input_name = f'{args.input_path}{data_name}'
    output_name = f'{args.output_path}res.mat'
    cmd = get_run_cmd(args, data_name, label_name)
    print(f'Run Command:{cmd}')
    os.system(cmd)
    print('Run Command end!')
    if not os.path.exists(output_name):
        print('Not find target file\nCheck error')
        return
    e1_meter = AverageMeter()
    e2_meter = AverageMeter()
    g1_meter = AverageMeter()
    if args.task != 4:
        truth_mat = loadmat(input_name)
        truth_mat = truth_mat['data']
        output_res = loadmat(output_name)
        E1 = np.squeeze(output_res['E1'])
        G1 = np.squeeze(output_res['G1'])
        e1_thd = e1_dict[args.task]
        g1_thd = g1_dict[args.task]
        for i in range(truth_mat.shape[0]):
            e1 = E1[i]
            g1 = G1[i]
            check = e1 < e1_thd and g1 < g1_thd
            if not check:
                print(
                    f'Find result exceed threshold\nTask:{args.task}\tIndex:{i}\tE1:{e1}\tG1:{g1}\nCheck error')
                return
            e1_meter.update(e1)
            g1_meter.update(g1)
    else:
        truth_mat = loadmat(input_name)
        truth_mat = truth_mat['data']
        output_res = loadmat(output_name)
        E2 = np.squeeze(output_res['E2'])
        e2_thd = e2_dict[args.task]
        for i in range(truth_mat.shape[0]):
            e2 = E2[i]
            check = e2 < e2_thd
            if not check:
                print(
                    f'Find result exceed threshold\nTask:{args.task}\tIndex:{i}\tE2:{e2}\nCheck error')
                return
            e2_meter.update(e2)
    print(
        f'Simple check pass!\nTask:{args.task}\tE1:{e1_meter.avg}\tG1:{g1_meter.avg}\tE2:{e2_meter.avg}')


if __name__ == '__main__':
    args = get_args()
    prepare(args)
    if args.measure:
        measure(args)
    elif args.complete_check:
        complete_check(args)
    elif args.simple_check:
        simple_check(args)
    clean_tmpfile(args)
