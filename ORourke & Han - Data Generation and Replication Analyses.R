#O'Rourke & Han (2023) "Considering the Distributional Form of Zeroes When Calculating Mediation Effects with Zero-Inflated Count Outcomes"
# R MplusAutomation for data generation and replication analysis scripts


install.packages("MplusAutomation")
library(MplusAutomation)


# Set up file paths to switch easily between home/work - choose 1
# Home file path
homepath <- "C:/Users/horourke/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/"
filepath <- homepath
# Work file path
workpath <- "C:/Users/HO Wexford/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/"
filepath <- workpath

#TEST PATH
testpath <- "C:/Users/HO Wexford/Desktop/test/"
filepath <- testpath


# Factors and n values to iterate over
factors <- c("ZINB", "ZIP")
n_values <- c(100, 250, 500, 750, 1500)
reps <- 3
dist <- "pi"
yvarspec <- "y"


# Nested loop over factors and n values
for (f in factors) {
  # Set the dist value based on the factor
  if (f == "ZINB") {
    dist <- "nbi"
    yvarspec <- "y"
    matrows <- 16
  } else if (f == "ZIP") {
    dist <- "pi"
    yvarspec <- "!y"
    matrows <- 15
  }
  for (n in n_values) {
    # Create the path for the script
    scr_path <- paste0(filepath, f, "/N = ", n)
    
    #Create directory for each combo of conditions
    if (!dir.exists(scr_path)) {
      dir.create(scr_path, recursive = TRUE)
    }
  
###############################    
# Write data generation scripts
###############################    
    
    # Open the file for writing
    dgscr <- paste0(scr_path, "/", f, "_n", n, ".inp")
    writeLines(sprintf('title: Data Generation Script for ORourke & Han MC Simulation, %s n=%d;
montecarlo:			
	names = x m y;
	seed = 53487;
	nobs = %d;
	nreps = %d;
	generate = y(%s);
    COUNT = y(%s);
    REPSAVE = ALL;
	save = %s%d_*.dat;
ANALYSIS: 
estimator=ml;
integration=montecarlo;

model population:
        [m@0];
        m@1;
        [x@0];
        x@0.25;
  		y ON x*.01 m*0.14;
   		[y@3];
        %s@1;
        y#1 ON x*-.01 m*-0.14;
  		[y#1@0];
        m ON x*0.59;
MODEL:
        [m] (m);
        m*1;
   		[y*3] (int_c);
        %s*1;
  		[y#1*0] (int_z);
  		y ON x*.01 (cpc);
        y ON m*0.14 (bc);
        y#1 ON x*-.01 (cpz);
        y#1 ON m*-0.14 (bz);
        m ON x*0.59 (a);
OUTPUT: TECH9;', 
                f, n, n, reps, dist, dist, f, n, yvarspec, yvarspec),
               con = dgscr)
    
    #Simulate data
    runModels(scr_path)
    
    #Create folder and file names
    folder_path <- paste0(filepath, f, "/N = ", n)
    file_suffix_n <- paste0(f, "_n", n)
    file_suffix_list <- paste0(tolower(f), n, "_list")
    file_suffix_r <- paste0(tolower(f), n, "_r")
    
    # Move _list files to datagen folder before creating .inp files
    datagen <- paste0(filepath, "datagen")
    if (!dir.exists(datagen)) {
      dir.create(datagen, recursive = TRUE)
    }
    file.copy(paste0(folder_path, "/", file_suffix_list, ".dat"), datagen)
    file.remove(paste0(folder_path, "/", file_suffix_list, ".dat"))
    
    file.copy(paste0(folder_path, "/", file_suffix_n, ".inp"), datagen)
    file.remove(paste0(folder_path, "/", file_suffix_n, ".inp"))
    file.copy(paste0(folder_path, "/", file_suffix_n, ".out"), datagen)
    file.remove(paste0(folder_path, "/", file_suffix_n, ".out"))
    
    ## Write MplusAutomation scripts ##
    # Open the file for writing
    mpf <- paste0(scr_path, "/", f, n, "_r.inp")
    
    # Write the Mplus script to the file
    writeLines(sprintf('[[init]]
iterators = sample;
sample = 1:%d;
filename = "%s_n%d_[[sample]].inp";
outputDirectory = "%s";
[[/init]]

TITLE: CIE calculations for %s n=%d rep [[sample]];
DATA:
    File is %s%d_[[sample]].dat;

  VARIABLE:
      NAMES = y m x;
      USEVARIABLES = y m x;
      count = y(%s);
ANALYSIS:
  estimator=ml;
  integration=montecarlo;
  !bootstrap=500;
MODEL:
          [m] (m);
          m*1;
            [y*3] (int_c);
          %s*1;
            [y#1*0] (int_z);
            y ON x*.01 (cpc);
          y ON m*0.14 (bc);
          y#1 ON x*-.01 (cpz);
          y#1 ON m*-0.14 (bz);
          m ON x*0.59 (a);

  MODEL CONSTRAINT:
  NEW(x_0_c x_1_c x_0_z_old x_0_z_li x_1_z_old x_1_z_li);

  x_0_c = (a*bc)*exp(int_c + bc*m + cpc*0);
  x_1_c = (a*bc)*exp(int_c + bc*m + cpc*1);

  x_0_z_old = a*bz*exp(int_z + bz*m + cpz*0);
  x_0_z_li = a*bz*(exp(int_z + bz*m + cpz*0)/((1+exp(int_z + bz*m + cpz*0))^2));

  x_1_z_old = a*bz*exp(int_z + bz*m + cpz*1);
  x_1_z_li = a*bz*(exp(int_z + bz*m + cpz*1)/((1+exp(int_z + bz*m + cpz*1))^2));
  OUTPUT: CINTERVAL(bcbootstrap);
  OUTPUT: TECH1;
  SAVEDATA: RESULTS ARE %s%d_est_[[sample]].dat;', 
                reps, f, n, scr_path, f, n, f, n, dist, yvarspec, f, n),
               con = mpf)

    # Create .inp files for all conditions and replications
    createModels(paste0(folder_path, "/", file_suffix_r, ".inp"))
    
    #Move createmodels .inp/.out files out of condition folder before running all models
    createm <- paste0(filepath,"/createmodels")
    if (!dir.exists(createm)) {
      dir.create(createm, recursive = TRUE)
    }
    file.copy(paste0(folder_path, "/", file_suffix_r, ".inp"), createm)
    file.remove(paste0(folder_path, "/", file_suffix_r, ".inp")) 

  #Run models for all reps
  runModels(folder_path)
  
  #############################
  # SET UP RESULTS FOR ANALYSIS
  #############################
    
    #Read in parameter estimates from all output files
    params <- paste0(f, "_", n, "_est")
    assign(params, readModels(
      target = folder_path,  
      recursive = TRUE,
      what = "parameters",
      quiet = TRUE
    ))
    
    #Initialize vector for reps
    estsname <- paste0(f, "_", n, "_ests_")
    ests_ <- assign(estsname, vector("list", reps))

    #Create empty matrix with estimates for all replications
    mat <- matrix(0, nrow = reps, ncol = matrows)
     
    for (i in 1:reps) {
      #Change filename depending on home/work location
      #filename <- paste0("C..Users.HO.Wexford.Dropbox..ASU..Papers.O.Rourke...MZI2.R.R.Analyses.", f, ".N...", n, ".", tolower(f), "_n", n, "_", i, ".out")
      filename <- paste0("C..Users.HO.Wexford.Desktop.test.", f, ".N...", n, ".", tolower(f), "_n", n, "_", i, ".out")
      ests_ <- get(params)[[filename]][["parameters"]][["unstandardized"]][["est"]]
      vec <- unlist(ests_)
      vec2 <- t(vec)
      mat[i, ] <- as.numeric(vec2)
    }
    #Create data frame from last list
    #ests <- paste0(f, "_", n, "_ests")
    #assign(ests, data.frame(t(sapply(ests_,c))))
    ests <- as.data.frame(mat)
    
    #Read out data as file for each condition to easily analyze later
    write.table(ests, file=paste0(folder_path, "/", file_suffix_n, "_ests.dat"), row.names=FALSE, sep="\t", quote=FALSE)
  }
}

#Create one dataset with all sim conditions for analysis
    length <- (length(factors)*length(n_values))
    simlist <- list()
    
    for (f in factors) {
      for (n in n_values) {
        simcond <- paste0(f,"_", n, "_ests")
        simlist[[simcond]] = get(simcond)
        }
      }
simests <- do.call(rbind,simlist)

#Create id, factor, sample size variables

# Initialize a row counter
row1 <- 1

# Loop over the conditions
for (f in factors) {
  for (n in n_values) {
    # Loop over the replications
    for (rep in 1:reps) {
      # Add a row to the dataframe
      simests[row1, "rep"] <- rep
      simests[row1, "factor"] <- f
      simests[row1, "n"] <- n
      
      # Increment the row counter
      row1 <- row1 + 1
    }
  }
}

namestr <- paste0(ZINB_100_est[["C..Users.horourke.Dropbox..ASU..Papers.O.Rourke...MZI2.R.R.Analyses.ZINB.N...100.ZINB_n100.out"]][["parameters"]][["unstandardized"]][["paramHeader"]],
                  " ",
                  ZINB_100_est[["C..Users.horourke.Dropbox..ASU..Papers.O.Rourke...MZI2.R.R.Analyses.ZINB.N...100.ZINB_n100.out"]][["parameters"]][["unstandardized"]][["param"]] 
                   )
print(namestr)
names(simests) <- c("cp_counts", "b_counts", "cp_zeroes", "b_zeroes", "a", "m_int", "y_z__int", "y_c_int", "m_resvar",
                    "dispers", "x_0_c", "x_1_c", "erp0", "erp1", "li0", "x_0_z_old", "x_0_z_li", "x_1_z_old",
                    "x_1_z_li", "rep", "factor", "n"
                    )
simests <- simests[,c(20,21,22,1:19)]
simests <- subset(simests, select = -c(erp0,erp1,li0))

#Save estimates data as one large .dat file
write.table(simests, file="C:/Users/horourke/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/all_ests.dat", row.names=FALSE, sep="\t", quote=FALSE)

attach(simests)

#Print tables of CIE means by sim condition
aggregate(cbind(simests$x_0_c, simests$x_1_c,
                simests$x_0_z_old, simests$x_0_z_li, 
                simests$x_1_z_old, simests$x_1_z_li), 
          list(simests$n, simests$factor), FUN=mean
          ) 

#Specify population values
nrows <- nrow(simests)
a_pop <- rep(0.5, nrows)
b_z_pop <- rep(-0.2, nrows)
b_c_pop <- rep(0.2, nrows)
cp_z_pop <- rep(-1.4, nrows)
cp_c_pop <- rep(1.4, nrows)
m_int_pop <- rep(0, nrows)
m_var_pop <- rep(1, nrows)
y_intz_pop <- rep(0, nrows)
y_intc_pop <- rep(3, nrows)
disp <- rep(1, nrows)

a <- simests$a
zero_b <- simests$b_zeroes
count_b <- simests$b_counts
count_cp <- simests$cp_counts
zero_cp <- simests$cp_zeroes
count_i2 <- simests$y_c_int
zero_i2 <- simests$y_z__int
m2_mean <- simests$m_int
med_c_0 = a * count_b * exp(count_i2 + (count_b * m2_mean) + (count_cp * 0))
med_c_1 = a * count_b * exp(count_i2 + (count_b * m2_mean) + (count_cp * 1))
exp0 <- exp(-(zero_i2 + (zero_b * m2_mean) + (zero_cp * 0)))
exp1 <- exp(-(zero_i2 + (zero_b * m2_mean) + (zero_cp * 1)))
new_med_z_0 = a * zero_b * (exp0 / ((1 + exp0)*(1 + exp0)))
new_med_z_1 = a * zero_b * (exp1 / ((1 + exp1)*(1 + exp1)))
old_med_z_0 = a * zero_b * exp(zero_i2 + (zero_b * m2_mean) + (zero_cp * 0))
old_med_z_1 = a * zero_b * exp(zero_i2 + (zero_b * m2_mean) + (zero_cp * 1))

#################
# Population CIEs
#################

#for  counts, x = 0:
med_c_0_pop = a_pop * b_c_pop * exp(y_intc_pop + (b_c_pop * m_int_pop) + (cp_c_pop * 0))

#for  counts, x = 1:
med_c_1_pop = a_pop * b_c_pop * exp(y_intc_pop + (b_c_pop * m_int_pop) + (cp_c_pop * 1))

#exponentials for zeroes
exp0_pop <- exp(-(y_intz_pop + (b_z_pop * m_int_pop) + (cp_z_pop * 0)))
exp1_pop <- exp(-(y_intz_pop + (b_z_pop * m_int_pop) + (cp_z_pop * 1)))

#for  zeroes, x = 0, old formula:
old_med_z_0_pop = a_pop * b_z_pop * exp0_pop

#for  zeroes, x = 1, old formula:
old_med_z_1_pop = a_pop * b_z_pop * exp1_pop

#for  zeroes, x = 0, new formula:
new_med_z_0_pop = a_pop * b_z_pop * (exp0_pop / ((1 + exp0_pop)*(1 + exp0_pop)))

#for  zeroes, x = 1, new formula:
new_med_z_1_pop = a_pop * b_z_pop * (exp1_pop / ((1 + exp1_pop)*(1 + exp1_pop)))

###pick up here calculating bias

#Calculate bias - edit to calculate  bias and population values using R formulas from checks in Dec
simests$a_bias_raw <- simests$a - a_pop
simests$a_bias_rel <- simests$a_bias_raw/a_pop
simests$b_z_bias_raw <- simests$b_zeroes - b_z_pop
simests$b_z_bias_rel <- simests$b_z_bias_raw/b_z_pop

simests$cie_z_0_bias_raw <- old_med_z_0 - new_med_z_0_pop
simests$cie_z_0_bias_rel <- simests$cie_z_0_bias_raw/new_med_z_0_pop

simests$cie_z_0_bias_raw2 <- new_med_z_0 - new_med_z_0_pop
simests$cie_z_0_bias_rel2 <- simests$cie_z_0_bias_raw2/new_med_z_0_pop

#Print tables of bias by sim condition
aggregate(cbind(simests$cie_z_0_bias_raw, simests$cie_z_0_bias_rel), list(simests$n, simests$factor), FUN=mean) 

aggregate(cbind(simests$cie_z_0_bias_raw2, simests$cie_z_0_bias_rel2), list(simests$n, simests$factor), FUN=mean) 

#Print tables of CIE means by sim condition
aggregate(cbind(simests$x_0_z_old, old_med_z_0_pop), 
          list(simests$n, simests$factor), FUN=mean
) 























##delete if not needed
# Move _list file to datagen folder before creating .inp files
file.copy(paste0(folder_path, file_suffix_list, ".dat"), "C:/Users/HO Wexford/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/datagen")
file.remove(paste0(folder_path, file_suffix_list, ".dat"))

# Create .inp files for all conditions and replications
createModels(paste0(folder_path, file_suffix_r, ".inp"))

#Move datagen and createmodels .inp/.out files out of condition folder
file.copy(paste0(folder_path, file_suffix_r, ".inp"), "C:/Users/HO Wexford/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/createmodels")
file.remove(paste0(folder_path, file_suffix_r, ".inp"))
file.copy(paste0(folder_path, file_suffix_n, ".inp"), "C:/Users/HO Wexford/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/datagen")
file.remove(paste0(folder_path, file_suffix_n, ".inp"))
file.copy(paste0(folder_path, file_suffix_n, ".out"), "C:/Users/HO Wexford/Dropbox (ASU)/Papers/O\'Rourke - MZI2/R&R/Analyses/datagen")


simests$id <- 1:nrow(simests)
simests$factor <- vector(, nrow(simests))
simests$factor[simests$id <= length(n_values)*reps] <- "ZINB"
simests$factor[simests$id > length(n_values)*reps] <- "ZIP"
simests$n <- vector(, nrow(simests))
simests$n[simests$id <= reps] <- 100
simests$n[simests$id > reps & simests$id <= reps*2] <- 250
simests$n[simests$id > reps*2 & simests$id <= reps*3] <- 500
simests$n[simests$id > reps*3 & simests$id <= reps*4] <- 750
simests$n[simests$id > reps*4 & simests$id <= reps*5] <- 1500
file.remove(paste0(folder_path, file_suffix_n, ".out"))


