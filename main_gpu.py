import torch
import torch.nn as nn
import torch.utils.data as data_utils
import models
from mmd import mix_rbf_mmd2
from utils import *
import math
import argparse


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
sigma_list = [1, 2, 4, 8, 16]


def train(original, batch_size, base_lr, lr_step, num_epochs, noise_level, lamda, num_iter):
    # network initialization
    dim_input = original.shape[1]
    net = models.network(dim_input)
    net.to(device)
    
    # make DataLoader
    original_tensor = torch.from_numpy(original).float().to(device)
    index_tensor = torch.arange(0, original.shape[0]).to(device)
    dataset_original = data_utils.TensorDataset(original_tensor, index_tensor)
    loader_original = data_utils.DataLoader(dataset_original, batch_size, shuffle=False)
    background = generate_bg(original, noise_level)
    background_tensor = torch.from_numpy(background).float().to(device)
    dataset_background = data_utils.TensorDataset(background_tensor, index_tensor)
    loader_background = data_utils.DataLoader(dataset_background, batch_size, shuffle=False)
    
    # start training
    print("start training...")
    for epoch in range(num_epochs):            
        learning_rate = base_lr / math.pow(2, math.floor(epoch / lr_step))
        recon_loss, dist_loss, loss = train_epoch(epoch, net, loader_original, loader_background, learning_rate, lamda)
        print("epoch ", epoch, " recon_loss: ", recon_loss.detach().numpy(), " dist_loss: ", dist_loss.detach().numpy(), 
              " total: ", loss.detach().numpy())
    return net
    
    
def train_epoch(epoch, net, data_loader1, data_loader2, learning_rate, lamda):
    
    optimizer = torch.optim.Adam([
        {'params': net.encoder1.parameters()},
        {'params': net.encoder2.parameters()},
        {'params': net.decoder1.parameters()},
        {'params': net.decoder2.parameters()}
    ], lr=learning_rate, weight_decay=1e-4)
    mse = nn.MSELoss()
    
    net.train()
    
    total_recon = torch.FloatTensor([0])
    total_dist = torch.FloatTensor([0])
    iterator1 = iter(data_loader1)
    iterator2 = iter(data_loader2)
    for i in range(len(data_loader1)):
        original, = next(iterator1)[:1]
        background, = next(iterator2)[:1]
        original_latent2, original_recon, _ , _ = net(original)
        background_latent2, _ ,decoded_background, _ = net(background)
        
        loss_recon = mse(original_recon, original) + mse(decoded_background, background)
        loss_dist = mix_rbf_mmd2(original_latent2, background_latent2, sigma_list)
        loss = loss_recon + lamda * loss_dist
        total_recon += loss_recon.cpu()
        total_dist += loss_dist.cpu()
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    
    avg_recon = total_recon / len(data_loader1)
    avg_dist = total_dist / len(data_loader1)
    avg_loss = avg_recon + avg_dist
    return avg_recon, avg_dist, avg_loss


def test(net, data):
    net.eval()
    _, _, _, denoised_patterns = net(torch.from_numpy(data).float().to(device))
    # ground_truth = torch.from_numpy(ground_truth).float().to(device)
    # diff = denoised_patterns - ground_truth
    # diff_nums = torch.sum(torch.abs(diff))
    # print("Wrong numbers in the denoised pattern matrix: ", diff_nums.cpu().detach().numpy())
    return denoised_patterns


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Training the DRLB model")
    parser.add_argument('-dataname', type=str, default=[], help="filename of your input matrix")
    parser.add_argument('-batch_size', type=int, default=8, help="batch sze for training")
    parser.add_argument('-base_lr', type=float, default=0.001, help="initial value of learning rate")
    parser.add_argument('-lr_step', type=float, default=10, help="step size for learning rate decay")
    parser.add_argument('-num_epochs', type=int, default=100, help="number of epochs in each iteration")
    parser.add_argument('-alpha', type=float, default=3, help="alpha value for background generation")
    parser.add_argument('-lamda', type=float, default=0.3, help="weight of the distribution loss")
    parser.add_argument('-num_iter', type=int, default=0, help="number of total iterations")

    args = parser.parse_args()
    dataname = args.dataname
    # dataname = 'scRNA-seq/GSE140228/top2k.csv'
    # dataname = 'simu7.5/data.csv'
    # gtname = 'simu4/patterns.csv'
    batch_size = args.batch_size
    base_lr = args.base_lr
    lr_step = args.lr_step
    num_epochs = args.num_epochs
    noise_level = args.alpha
    lamda = args.lamda
    num_iter = args.num_iter
    
    binary_matrix = read_data(dataname).astype(float)
    # binary_matrix = np.transpose(binary_matrix)
    # binary_matrix = binarize(binary_matrix)
    # ground_truth = read_data(gtname)
    # ground_truth = binary_matrix
    network = train(binary_matrix, batch_size, base_lr, lr_step, num_epochs, noise_level, lamda, num_iter)    
    denoised_pattern_matrix = test(network, binary_matrix)
    denoised_pattern_matrix = denoised_pattern_matrix.cpu().detach().numpy()
    # denoised_pattern_matrix = np.transpose(denoised_pattern_matrix)
    
    df = pd.DataFrame(denoised_pattern_matrix)
    df.to_csv('debiased.csv')
    # df.to_csv('scRNA-seq/GSE140228/denoised2.csv')
