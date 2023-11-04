#O'Rourke & Han (2023) "Considering the Distributional Form of Zeroes When Calculating Mediation Effects with Zero-Inflated Count Outcomes"
# R MplusAutomation for data generation and replication analysis scripts
# Ensure you have administrative permission to copy and remove files before running


install.packages("MplusAutomation")
library(MplusAutomation)

# Set up file paths to switch easily between home/work - choose 1
# "filepath" will be the directory to the folder where all simulated data will be stored

# Home file path
homepath <- "C:/myfiles/home/"
filepath <- homepath
# Work file path
workpath <- "C:/myfiles/work"
filepath <- workpath


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
    
###############################    
# Write MplusAutomation scripts
###############################    
    
    # Open the file for writing
    mpf <- paste0(scr_path, "/", f, n, "_r.inp")
    
    # Write the Mplus script to the file
    writeLines(sprintf('[[init]]
iterators = sample;
sample = 1:%d;
filename = "%s_n%d_[[sample]].inp";
outputDirectory = "%s";
[[/init]]

TITLE: ORourke & Han MC Simulation, Replication Analyses for %s n=%d rep [[sample]];
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
      #Dot vs. slash convention is for MplusAutomation
      #filename <- paste0("C..myfiles.home.", f, ".N...", n, ".", tolower(f), "_n", n, "_", i, ".out")
      filename <- paste0("C..myfiles.work.", f, ".N...", n, ".", tolower(f), "_n", n, "_", i, ".out")
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

names(simests) <- c("cp_counts", "b_counts", "cp_zeroes", "b_zeroes", "a", "m_int", "y_z_int", "y_c_int", "m_resvar",
                    "dispers", "x_0_ab_c", "x_1_ab_c", "erp0", "erp1", "li0", "x_0_ab_ll", "x_0_ab_lg", "x_1_ab_ll",
                    "x_1_ab_lg", "rep", "factor", "n"
                    )
simests <- simests[,c(20,21,22,1:19)]
simests <- subset(simests, select = -c(erp0,erp1,li0))

#Save estimates data as one large .dat file
write.table(simests, file=paste0(filepath, "all_ests.dat"), row.names=FALSE, sep="\t", quote=FALSE)

