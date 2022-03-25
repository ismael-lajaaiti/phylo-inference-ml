
# Import libraries 
library(Hmisc)
library(corrplot)
library(ape)


# Load Predictions from file 
pred.list <- readRDS("params-testset-bisse.RDS")
pred.df <- pred.list %>% as.data.frame() # convert to data.frame

pred.df.lambda <- pred.df[c(1,3,5,7,9,11,13)]  # data.frame for lambda00
pred.df.q      <- pred.df[c(2,4,6,8,10,12,14)] # data.frame for q01

# Compute the prediction error by substracting the 1st column (true values)
# to the other columns (predicted values)
error.df.lambda <- pred.df.lambda[2:7] - pred.df.lambda[,1] 
error.df.q      <- pred.df.q[2:7]      - pred.df.q[,1]

# Rename columns
column_names <- c("MLE", "rnn-ltt", "cnn-ltt", "cnn-enc", "dnn-ss",
                  "gnn-phy")
colnames(error.df.lambda) <- column_names
colnames(error.df.q)      <- column_names

# Compute correlations matrices
cor.lambda <- cor(error.df.lambda, method = c('spearman'))
cor.q      <- cor(error.df.q,      method = c('spearman'))
rcorr.lambda <- rcorr(as.matrix(error.df.lambda), type = c('spearman'))$P %>% as.matrix()

par(mfrow=c(1,2))

p.mat.q <- cor.mtest(error.df.q, conf.level = 0.95)$p
p.mat.lambda <- cor.mtest(error.df.lambda, conf.level = 0.95)$p
corrplot(cor.lambda, p.mat = p.mat.lambda, sig.level = 0.05, order = 'hclust',
         addrect = 3, title = "lambda", tl.cex = .8, tl.pos = "tl",
         tl.col = "black", cl.pos ="n", tl.srt = 45, tl.offset = .5,
         mar = c(5, 4, 5, 1))


corrplot(cor.q, p.mat = p.mat.q, sig.level = 0.05, order = 'hclust',
         addrect = 2, title = "q", tl.cex = .8, tl.pos = "t", 
         tl.col = "black", mar = c(5, 3, 5, 2), tl.srt = 45)

corrplot.mixed(cor.lambda, order = 'hclust', addrect = 3)

corrplot(cor.lambda, method = c("circle"), type = c("full"), order = c("hclust"), tl.cex = .8)
corrplot(cor.q, method = c("circle"), type = c("full"), order = c("hclust"), tl.cex = .8)

palette = colorRampPalette(c("white", "green")) (20)
heatmap(x = cor.lambda, col = palette, symm = TRUE)


# Empirical phylogenies 
data("bird.families")  # 137 tips
data("chiroptera")     # 916 tips
data("Phyllostomidae") # 150 tips
data("Cetacea")        # 87  tips 


Cetacea
chiroptera
Phyllostomidae
bird.families
plot(Phyllostomidae)
