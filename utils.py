import pandas as pd
import numpy as np
import torch


def read_data(file):
    df = pd.read_csv(file)
    data = df.values[:, 1:]
    return data


def binarize(data):
    return (data>0)*1


# generate a background matrix with a specified noise level
def generate_bg(data, noise_level):
    n_rows = data.shape[0]
    n_cols = data.shape[1]
    rows_prob = np.zeros((n_rows))
    cols_prob = np.zeros((n_cols))
    for i in range(n_rows):
        rows_prob[i] = sum(data[i, :]) / n_cols
    for j in range(n_cols):
        cols_prob[j] = sum(data[:, j]) / n_rows
    bg = np.zeros_like(data)
    for i in range(n_rows):
        for j in range(n_cols):
            prob = rows_prob[i] * cols_prob[j] * noise_level
            if prob > 1:
                prob = 1  # added in case prob > 1
            bg[i, j] = np.random.binomial(1, prob)
    return bg


def generate_bg2(data):
    n_rows = data.shape[0]
    n_cols = data.shape[1]
    rows_prob = np.zeros((n_rows))
    cols_prob = np.zeros((n_cols))
    for i in range(n_rows):
        rows_prob[i] = sum(data[i]) / n_cols
    for j in range(n_cols):
        cols_prob[j] = sum(data[:, j]) / n_rows
    bg = np.zeros_like(data)
    for i in range(n_rows):
        for j in range(n_cols):
            row_hit = np.random.binomial(1, rows_prob[i])
            col_hit = np.random.binomial(1, cols_prob[j])
            bg[i, j] = ((row_hit + col_hit) > 0)*1
    return bg


# convert a file with format "user | item | rate | timestamp" to 0-1 matrix
def udata2mat(path):
    users = []
    items = []
    rates = []
    file = open(path)
    lines = file.readlines()
    for line in lines:
        line = line.split('\t')
        users.append(int(line[0]))
        items.append(int(line[1]))
        rates.append(float(line[2]))
    n_users = max(users)
    n_items = max(items)
    mat = np.zeros((n_users, n_items))
    for i in range(len(users)):
        if rates[i] > 3:
            mat[users[i]-1][items[i]-1] = 1
    return mat


# update background based on estimated patterns
def update_bg(data, patterns):
    mask = 1 - patterns
    bg = torch.mul(data, mask)
    return bg.cpu().detach().numpy()