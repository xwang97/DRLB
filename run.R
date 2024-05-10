# 1. load the dataset
pat_num = 2
data <- as.matrix(read.csv("data/simu1/data.csv")[, 2:501])
patterns <- as.matrix(read.csv("data/simu1/patterns.csv")[, 2:501])
background <- as.matrix(read.csv("data/simu1/bg.csv")[, 2:501])

# 2. run BMF without debiasing and evaluate
result <- MEBF(data, DIM = pat_num, Thres = 0.6)
# result <- ASSO(data, DIM = pat_num, Thres = 0.6)
# result <- PANDA(data, DIM = pat_num)
recon_error <- pat_recon_err(result[[1]], result[[2]], patterns)
print(paste0("without debiasing: ", recon_error[[1]], "   ", recon_error[[2]]))

# 3. debias with BIND
mat_after_bind <- denoise_with_bind(data)
myplot(mat_after_bind, patterns)
result <- MEBF(mat_after_bind, DIM = pat_num, Thres = 0.6)
# result <- ASSO(mat_after_bind, DIM = pat_num, Thres = 0.6)
# result <- PANDA(mat_after_bind, DIM = pat_num)
recon_error <- pat_recon_err(result[[1]], result[[2]], patterns)
print(paste0("after BIND: ", recon_error[[1]], "   ",recon_error[[2]]))

# 4. evaluate the result of DRLB
mat_after_nn <- as.matrix(read.csv("debiased.csv")[, 2:501])
myplot(mat_after_nn, patterns)
result <- MEBF(mat_after_nn, DIM = pat_num, Thres = 0.6)
# result <- ASSO(mat_after_nn, DIM = pat_num, Thres = 0.6)
# result <- PANDA(mat_after_nn, DIM = pat_num)
recon_error <- pat_recon_err(result[[1]], result[[2]], patterns)
print(paste0("after DRLB: ", recon_error[[1]], "   ",recon_error[[2]]))


