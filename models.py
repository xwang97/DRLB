import torch
import torch.nn as nn
from torch.autograd import Function


def init_weights(m):
    """ initialize weights of fully connected layer
    """
    if type(m) == nn.Linear:
        nn.init.xavier_uniform_(m.weight)
        m.bias.data.fill_(0.001)


def init_uniform(m):
    """
    initialize weights with uniform distributions
    """
    if type(m) == nn.Linear:
        nn.init.uniform_(m.weight, 0, 0.2)
        nn.init.uniform_(m.bias, 0, 0.01)


class round_func(Function):
    @staticmethod
    def forward(ctx, input):
        ctx.input = input
        return (input>0.6).float()

    @staticmethod
    def backward(ctx, grad_output):
        grad_input = grad_output.clone()
        return grad_input

class encoder(nn.Module):
    def __init__(self, dim_inputs):
        super(encoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(dim_inputs, 200),
            nn.ReLU(),
            nn.Linear(200, 20)
        )
        self.encoder.apply(init_uniform)
    
    def forward(self, x):
        x = self.encoder(x)
        return x


class decoder(nn.Module):
    def __init__(self, dim_inputs):
        super(decoder, self).__init__()
        self.decoder = nn.Sequential(
            nn.Linear(20, 200),
            nn.ReLU(),
            nn.Linear(200, dim_inputs),
            #nn.Sigmoid()
            nn.ReLU()
        )
        self.decoder.apply(init_uniform)
    
    def forward(self, x):
        x = self.decoder(x)
        return x


class network(nn.Module):
    def __init__(self, dim_inputs):
        super(network, self).__init__()
        self.encoder1 = encoder(dim_inputs)
        self.encoder2 = encoder(dim_inputs)
        self.decoder1 = decoder(dim_inputs)
        self.decoder2 = decoder(dim_inputs)
    
    def forward(self, x):
        latent_patterns = self.encoder1(x)
        latent_background = self.encoder2(x)
        pattern_matrix = round_func.apply(self.decoder1(latent_patterns))
        background_matrix = round_func.apply(self.decoder2(latent_background))
        recon_matrix = round_func.apply(pattern_matrix + background_matrix)
        return latent_background, recon_matrix, background_matrix, pattern_matrix
        
