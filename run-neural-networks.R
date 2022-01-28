#'------------------------------------------------------------------------------
#'------------------------RUNNING NEURAL NETWORKS-------------------------------
#'------------------------------------------------------------------------------




run_dnn_summary_statistics <- function(trees_train, trees_valid, trees_test,
                                       param.range, n_taxa, epoch_max = 100,
                                       n_hidden = 100, n_layer = 4,
                                       load_ss = TRUE){
  
  
  
  
  
}
  
  
  
  
  
  error ~ true * model 


fit <- lm(Ozone ~ Temp, data = airquality)
library(effects)
plot(allEffects(fit))


library(JuliaCall)
julia <- julia_setup()
julia_library("PANDA")

x <- 2 
julia_eval(paste("sqrt(", x, ")", sep = ""))

julia_command("tree = sim_ClaDS2_ntips(100,0.1,1.,0.5,1.)")
julia_command("infer_ClaDS(tree, 1, 1.05, Inf, 1., 0.25, 10, 50, 100)")


library(RPANDA)

obj = sim_ClaDS(lambda_0=1.,    
                mu_0=0.0,      
                sigma_lamb=0.1,         
                alpha_lamb=1.1,     
                condition="taxa",    
                taxa_stop = 100,    
                prune_extinct = TRUE)  

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

plot_ClaDS_phylo(tree,speciation_rates)

sampler = fit_ClaDS0(tree=tree,
                     name=NULL,
                     iteration=500000,
                     thin=2000,
                     update=1000, adaptation=5) 




MAPS = getMAPS_ClaDS0(tree, sampler, burn = 1/2, thin = 10)
MAPS[1:3]
