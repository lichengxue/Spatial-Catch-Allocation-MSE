# devtools::install_github("lichengxue/wham", dependencies=TRUE, ref="lab", lib = "~/R/x86_64-redhat-linux-gnu-4.3.2/")

setwd("/gpfs/scratch/chengxli/MBI")
.libPaths("~/R/x86_64-redhat-linux-gnu-4.3.2/")

nsim <- 1
models_to_run <- paste0("Mod", 1)  # "Mod1"–"Mod6"

library(wham)
library(whamMSE)

F_est <- readRDS("F_est.RDS")
waa <- readRDS("WAA.RDS")

n_ages    <- 8
years     <- 1989:2023
MSE_years <- 50
n_seasons <- 11
n_stocks  <- 2
n_regions <- 2
n_fleets  <- 4
n_indices <- 2
Fbar_ages <- 8
fracyr_seasons <- c(rep(1,5),2,rep(1,5))/12

user_waa <- list()
user_waa$waa <- array(NA, dim = c(8,length(years)+MSE_years,8))
for (age in 1:8) {
  for (indx in 1:8) {
    user_waa$waa[indx,,age] = rep(waa[indx,age], length(years)+MSE_years)
  }
}
user_waa$waa_pointer_fleets  <- 1:4
user_waa$waa_pointer_indices <- 5:6
user_waa$waa_pointer_ssb     <- 7:8
user_waa$waa_pointer_M       <- 7:8

user_maa <- rbind(c(0,0.51,0.95,1,1,1,1,1),
                  c(0,0.45,0.85,0.95,0.99,1,1,1))

q             <- c(7.370412e-05, 2.832577e-05)
fleet_regions <- c(1,1,2,2)
catch_Neff    <- c(50,50,50,50)
catch_cv      <- c(0.05,0.15,0.05,0.15)
index_regions <- c(1,2)
index_Neff    <- 50
index_cv      <- 0.4

info <- generate_basic_info(
  base.years       = years, 
  n_feedback_years = MSE_years,
  n_stocks         = n_stocks,
  n_regions        = n_regions,
  n_indices        = n_indices,
  n_fleets         = n_fleets,
  n_ages           = n_ages,
  n_seasons        = n_seasons,
  catch_info       = list(catch_cv = catch_cv, catch_Neff = catch_Neff),
  index_info       = list(index_cv = index_cv, index_Neff = index_Neff, q = q, fracyr_indices = 0.25),
  fracyr_seasons   = fracyr_seasons,
  fracyr_spawn     = 0.5,
  user_waa         = user_waa,
  user_maturity    = user_maa,
  move_dyn = 1
)
basic_info <- info$basic_info
catch_info <- info$catch_info
index_info <- info$index_info
F_info     <- info$F

F_info$F[1:35,] <- F_est

basic_info$NAA_where <- array(1, dim = c(n_stocks,n_regions,n_ages))
basic_info$NAA_where[1,2,1] <- 0
basic_info$NAA_where[2,1,]  <- 0

move <- list(stock_move = c(TRUE,FALSE), separable = TRUE)
move$must_move <- array(0,dim = c(2,11,2))

move$can_move  <- array(0, dim = c(2,11,2,2))
move$can_move[1, 1:5 , 2, 1]  <- 1
move$can_move[1, 7:11, 1, 2]  <- 1

mean_vals <- array(0, dim = c(2,length(fracyr_seasons),2,1))
mean_vals[1,1:11,1,1] <- 0.02214863
mean_vals[1,1:11,2,1] <- 0.3130358

move$mean_vals  <- mean_vals
move$mean_model <- matrix("stock_constant", 2,1)

move$year_re <- matrix("iid",2,1)

move$sigma_vals = 0.5
move$prior = array(0,dim= c(2,11,2,1))

sigma_vals   <-  array(0.2, dim = c(n_stocks, n_regions, n_ages))
sigma_vals[,,1] = 0.75

NAA_re <- list(N1_model=rep("age-specific-fe",n_stocks),
               sigma=rep("rec+1",n_stocks),
               cor=rep("iid",n_stocks),
               recruit_model = 2,
               recruit_pars = list(16000,21000), 
               sigma_vals = sigma_vals)

M <- list(model = "constant", initial_means = array(0.4, dim = c(n_stocks, n_regions, n_ages)))

model <- c("age-specific","age-specific","logistic","logistic","age-specific","age-specific")
initial_pars_f1 <- c(0.031, 0.421, 0.858, 1, 1, 1, 1, 1)
initial_pars_f2 <- c(0.036, 0.355, 0.573, 0.782, 0.868, 0.928, 1, 1)
initial_pars_f3 <- c(2.381,0.402)
initial_pars_f4 <- c(2.788,0.849)
initial_pars_s1 <- c(0.091, 0.463, 0.905, 0.962, 1, 1, 1, 1)
initial_pars_s2 <- c(0.362, 1, 1, 1, 1, 1, 1, 1)

sel <- list(model=model,
            initial_pars=c(list(initial_pars_f1), list(initial_pars_f2), 
                           list(initial_pars_f3), list(initial_pars_f4),
                           list(initial_pars_s1), list(initial_pars_s2)),
            fix_pars = list(4:8, 7:8, NULL, NULL, 5:8, 2:8))

input <- prepare_wham_input(basic_info = basic_info, 
                            NAA_re = NAA_re,
                            move = move, 
                            M = M, 
                            selectivity = sel,
                            catch_info = catch_info, 
                            index_info = index_info, 
                            F = F_info
)

waa_info <- info$par_inputs$user_waa
input <- update_waa(input, waa_info = waa_info)

ini.NAA <- matrix(NA, n_ages, n_stocks)
ini.NAA[,1] <- c(4531.15069,2876.48828,1484.19315,692.67290,314.73335,142.50149,64.41915,52.58156)
ini.NAA[,2] <- c(16023.94518,10167.14058,5683.71704,2523.02265,972.93721,353.59712,125.80803,68.66246)

input$par$log_N1[] = 0
for (i in 1:n_regions) {
  input$par$log_N1[i,i,] = log(ini.NAA[,i])
}

random = input$random # check what processes are random effects
input$random = NULL # so inner optimization won't change simulated RE
om <- fit_wham(input, do.fit = F, do.brps = F, MakeADFun.silent = TRUE)

saveRDS(om,    "OM_MBI.RDS")
saveRDS(input, "input.RDS")
saveRDS(info,  "info.RDS")

## ==========================================
## 2. Parallel MSE Runs for Mod1–Mod6
## ==========================================

library(doParallel)
library(foreach)
options(connectionObserver = NULL)

n_cores <- min(length(nsim), max(1L, parallel::detectCores() - 1L))
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

# Setup cluster environment
parallel::clusterCall(cl, function() {
  .libPaths(c("~/R/x86_64-redhat-linux-gnu-4.3.2", .libPaths()))
  Sys.setenv(OMP_NUM_THREADS = 1,
             OPENBLAS_NUM_THREADS = 1,
             MKL_NUM_THREADS = 1,
             OMP_THREAD_LIMIT = 1)
  setwd("/gpfs/scratch/chengxli/MBI")
  invisible(NULL)
})

parallel::clusterEvalQ(cl, {
  library(wham)
  library(whamMSE)
})

# Model configurations
model_configs <- list(
  Mod1 = list(catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL,      weight_years = 3)),
  Mod2 = list(catch_alloc = list(weight_type = 1, method = "gear_equal", user_weights = NULL, weight_years = 3)),
  Mod3 = list(catch_alloc = list(weight_type = 2, method = "fleet_catch", user_weights = NULL, weight_years = 3)),
  Mod4 = list(catch_alloc = list(weight_type = 3, method = "index_equal", user_weights = NULL, weight_years = 3)),
  Mod5 = list(catch_alloc = list(weight_type = 3, method = "index_gear", user_weights = NULL,  weight_years = 3)),
  Mod6 = list(catch_alloc = list(weight_type = 4, method = "user_defined_regions", user_weights = c(0.1,0.9), weight_years = 3))
)

if (!dir.exists("Results")) dir.create("Results", recursive = TRUE)

n_ages <- 8

# Outer loop: loop over models Mod1–Mod6
for (mod_name in models_to_run) {
  cat("=== Running", mod_name, "===\n")
  
  foreach(i = nsim,
          .packages = c("wham","whamMSE"),
          #.export   = c("n_ages"),
          .inorder  = FALSE,
          .options.snow = list(preschedule = FALSE)) %dopar% {
            
            assign("n_ages", 8, envir = .GlobalEnv)
            
            om <- readRDS("OM_MBI.RDS")
            random <- c("log_NAA","mu_re")
            om_with_data <- update_om_fn(om, seed = 5201314 + as.integer(i), random = random)
            
            sel_em <- list(
              model=c("age-specific","age-specific","age-specific"),
              initial_pars=c(
                list(c(0.5,0.5,0.5,1,1,1,1,1)),
                list(c(0.5,0.5,0.5,0.5,1,1,1,1)),
                list(c(0.5,1,1,1,1,1,1,1))
              ),
              fix_pars = list(4:8, 5:8, 2:8)
            )
            
            sigma_vals <- array(0.2, dim = c(1, 1, 8))
            sigma_vals[,,1] <- 0.75
            
            NAA_re_em <- list(
              N1_model = "age-specific-fe",
              sigma="rec+1",
              cor="iid",
              recruit_model = 2,
              recruit_pars = list(37000),
              sigma_vals = sigma_vals
            )
            
            M_em <- list(model="constant", initial_means=array(0.4, dim = c(1,1,8)))
            
            assess.interval <- 2
            base.years <- 1989:2023
            terminal.year <- tail(base.years,1)
            assess.years <- seq(terminal.year, 2023+50-2, by = assess.interval)
            
            info <- readRDS("info.RDS")
            catch_alloc_config <- model_configs[[mod_name]]$catch_alloc
            
            n_ages = 8
            
            mod <- loop_through_fn(
              om = om_with_data,
              NAA_re_em = NAA_re_em,
              sel_em = sel_em,
              M_em = M_em,
              move_em = NULL,
              em_info = info,
              em.opt = list(separate.em = TRUE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
              aggregate_catch_info = list(n_fleets = 2, fleet_pointer = c(1,2,1,2), use_catch_weighted_waa = TRUE, 
                                          catch_Neff = c(50,50), catch_cv = c(0.05,0.15)),
              aggregate_index_info = list(n_indices = 1, index_pointer = c(1,1), use_catch_weighted_waa = TRUE, 
                                          index_Neff = 50, index_cv = 0.4),
              aggregate_weights_info = list(
                ssb_waa_weights = list(index = TRUE, fleet = FALSE, pointer = c(1,2)),
                maturity_weights = list(index = TRUE, fleet = FALSE, pointer = c(1,2))
              ),
              assess_years = assess.years,
              assess_interval = assess.interval,
              base_years = base.years,
              year.use = 35,
              add.years = TRUE,
              catch_alloc = catch_alloc_config,
              save.last.em = FALSE,
              seed = 5201314 + as.integer(i),
              FXSPR_init = 0.5
            )
            
            saveRDS(mod, file.path("Results", sprintf("%s_%03d.RDS", mod_name, as.integer(i))))
            NULL
          }
}

stopCluster(cl)
