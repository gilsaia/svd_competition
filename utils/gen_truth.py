import argparse
import scipy
import numpy as np
import os
from scipy.io import loadmat, savemat


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str, required=True)
    parser.add_argument('--input_name', type=str, required=True)
    parser.add_argument('--output_path', type=str, required=True)
    parser.add_argument('--svd', action='store_true')
    parser.add_argument('--inv', action='store_true')
    parser.add_argument('--transpose', action='store_true')
    return parser.parse_args()


def compute_transpose(args):
    input_path = f'{args.input_path}{args.input_name}'
    output_path = f'{args.output_path}{args.input_name}'
    input = loadmat(input_path)
    for key in input:
        if isinstance(input[key], np.ndarray):
            content = input[key]
            break
    if len(content.shape) == 3:
        changedim = (2, 0, 1)
    else:
        changedim = (1, 0)
    content = np.transpose(content, changedim)
    print(content[0, 0])
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    savemat(output_path, {'data': content})


def compute_svd(args):
    input_path = f'{args.input_path}{args.input_name}'
    output_path = f'{args.output_path}{args.input_name}'
    input = loadmat(input_path)
    for key in input:
        if isinstance(input[key], np.ndarray):
            content = input[key]
            break
    u, s, vh = np.linalg.svd(content)
    v = np.transpose(vh, (0, 2, 1))
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    savemat(output_path, {'u': u, 'v': v, 's': s})


if __name__ == '__main__':
    args = get_args()
    if args.transpose:
        compute_transpose(args)
    if args.svd:
        compute_svd(args)
