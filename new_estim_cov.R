# Created by Oleksandr Sorochynskyi
# On 27 08 18




source('~/IdV/R/conditional.R')

# The code to use for different stimuli types
# 1 = fullfield
# 2 = barmovie
# 3 = checkerboard

estim_cov <- function(data_stimulus, parameter_stimulus) {
  stopifnot(length(data_stimulus) == 1,
            length(parameter_stimulus) == 1,
            any(data_stimulus == c(1,2,3)),
            any(parameter_stimulus == c(1,2,3))
            )
  
  binner <- switch(data_stimulus,
                   bin_fullfield,
                   bin_barmovie,
                   bin_checkerboard) 
  data_stimulus <- switch(data_stimulus,
                         "fullfield",
                         "bar",
                         "check")
  parameter_stimulus <- switch(parameter_stimulus,
                         "fullfield",
                         "bar",
                         "check")
  
  filename <- paste("cov/",
                    data_stimulus,
                    "_data_",
                    parameter_stimulus,
                    "_param.Rdata", sep="");
  
  
  n0 <- 25
  n1 <- 19
  n <- n0+n1
  m0 <- binner(1:n0, cell.type= 1);
  m1 <- binner(1:n1, cell.type= 6);
  m <- array(NA, dim=c(n1+n0,dim(m0)[2],dim(m0)[3]))
  m[1:n0,,] <- m0;
  m[-(1:n0),,] <- m1;
  
  n.rep <- 1000;
  load("summary/params.RData");
  par <- list(model2.2=get(parameter_stimulus,params$model2.2$par),
              model3=get(parameter_stimulus,params$model3$par));
  
  noise_cov_real <- matrix(NA, ncol=n, nrow=n);
  noise_cov_simu <- matrix(NA, ncol=n, nrow=n);
  stimul_cov_real <- matrix(NA, ncol=n, nrow=n);
  stimul_cov_simu <- matrix(NA, ncol=n, nrow=n);
  
  pb <- txtProgressBar(min=0, n*(n-1), style= 3);
  
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      tryCatch({
        rand <- cond_model2.2_rand(i, j, m, par$model2.2[i,j], n.rep);
        
        noise_cov_real[i,j] <- cond_noise_cov(i,j, m);
        stimul_cov_real[i,j] <- cond_stimul_cov(i, j, m);
        
        noise_cov_simu[i,j] <- cond_noise_cov(1,2,rand);
        stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
      },
      error = function(e) {
        warning(paste("Could't fit model 2.2 for pair:",i,",",j,"."))
        
        noise_cov_real[i,j] <- NA;
        stimul_cov_real[i,j] <- NA;
        
        noise_cov_simu[i,j] <- NA;
        stimul_cov_simu[i,j] <- NA;
      });
      
      
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1);
    }
  }
  
  total_cov_real <- stimul_cov_real + noise_cov_real;
  total_cov_simu <- stimul_cov_simu + noise_cov_simu;
  
  model2.2_cov_real <- list(total=total_cov_real,
                          noise=noise_cov_real,
                          stimul=stimul_cov_real)
  model2.2_cov_simu <- list(total=total_cov_simu,
                          noise=noise_cov_simu,
                          stimul=stimul_cov_simu)
  
  # Now model3
  noise_cov_real <- matrix(NA, ncol=n, nrow=n);
  noise_cov_simu <- matrix(NA, ncol=n, nrow=n);
  stimul_cov_real <- matrix(NA, ncol=n, nrow=n);
  stimul_cov_simu <- matrix(NA, ncol=n, nrow=n);
  
  
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      tryCatch({
        rand <- cond_model3_rand(i, j, m, par$model3[i,j], n.rep);
        
        noise_cov_real[i,j] <- cond_noise_cov(i,j,m);
        stimul_cov_real[i,j] <- cond_stimul_cov(i, j, m);
        
        noise_cov_simu[i,j] <- cond_noise_cov(1,2,rand);
        stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
      },
      error = function(e) {
        warning(paste("Could't fit model 2.2 for pair:",i,",",j,"."))
        
        noise_cov_real[i,j] <- NA;
        stimul_cov_real[i,j] <- NA;
        
        noise_cov_simu[i,j] <- NA;
        stimul_cov_simu[i,j] <- NA;
      });
      
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1);
    }
  }
  
  close(pb);
  
  total_cov_real <- stimul_cov_real + noise_cov_real;
  total_cov_simu <- stimul_cov_simu + noise_cov_simu;
  
  model3_cov_real <- list(total=total_cov_real,
                          noise=noise_cov_real,
                          stimul=stimul_cov_real)
  model3_cov_simu <- list(total=total_cov_simu,
                          noise=noise_cov_simu,
                          stimul=stimul_cov_simu)
  
  model2.2 <- list(real=model2.2_cov_real,
                 simu=model2.2_cov_simu);
  model3 <- list(real=model3_cov_real,
                 simu=model3_cov_simu);
  
  var.name <- paste(data_stimulus, "_data_", parameter_stimulus, "_param_cov", sep="");
  assign(var.name, list(model2.2= model2.2,
                        model3= model3));
  
  saveit <- function(..., string, file) {
    x <- list(...)
    names(x) <- string
    save(list=names(x), file=file, envir=list2env(x))
  }
  
  saveit(var = get(var.name), string= var.name, file= filename)
  get(var.name)
}

# data
# 3|x x x
# 2|x x x
# 1|x_x_x
#   1|2|3 param

estim_cov(1,1); # done
estim_cov(1,2); # done
estim_cov(1,3); # done
estim_cov(2,1); # done
estim_cov(2,2); # done
estim_cov(2,3); # done
estim_cov(3,1); # done
estim_cov(3,2); # done
estim_cov(3,3); # done

# First load all the resulting files 
# Compile them together.
bar_data <- list(fullfield_param=bar_data_fullfield_param_cov,
                 bar_param=bar_data_bar_param_cov,
                 check_param=bar_data_check_param_cov);
check_data <- list(fullfield_param=check_data_fullfield_param_cov,
                 bar_param=check_data_bar_param_cov,
                 check_param=check_data_check_param_cov);
fullfield_data <- list(fullfield_param=fullfield_data_fullfield_param_cov,
                 bar_param=fullfield_data_bar_param_cov,
                 check_param=fullfield_data_check_param_cov);
covar <- list(bar_data=bar_data,
                check_data=check_data,
                fullfield_data=fullfield_data);

save(covar, file="summary/covar.RData");
