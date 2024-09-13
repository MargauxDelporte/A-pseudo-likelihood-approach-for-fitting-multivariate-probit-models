library(mvProbit)  # For mvProbit function
library(haven)     # For read_sas function
library(readxl)    # For read_excel function
library(xlsx)
setwd('')

################################
###create a function to run it #
################################
analyze_dataset=function(i){
  while(length(dir(tempdir()))>100){
    # Remove temporary files 
    tmp_dir <- tempdir()
    files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    file.remove(files)
  }
  # Initial parameters
  startv=c(rep(0,7),
           rep(0,7),
           rep(0,7),
           rep(0,3))
  parameters=c('intercept_1',
               'beta_sex_1',
               'beta_age_1',
               'beta_reg2_1',
               'beta_reg3_1',
               'beta_reg4_1',
               'beta_reg5_1',
               'intercept_2',
               'beta_sex_2',
               'beta_age_2',
               'beta_reg2_2',
               'beta_reg3_2',
               'beta_reg4_2',
               'beta_reg5_2',
               'intercept_3',
               'beta_sex_3',
               'beta_age_3',
               'beta_reg2_3',
               'beta_reg3_3',
               'beta_reg4_3',
               'beta_reg5_3',
               'r12',
               'r13',
               'r23')
  
  # Construct the file name dynamically
  file_path <- gsub(" ", "",paste("simulated_", i, ".xlsx"))
  
  # Read in the dataset
  dataset <- read_excel(file_path, NULL)
  
  
  #fit model
  formula=as.formula(cbind(y1b, y2b,y3b) ~ sex+age+as.factor(region))
  result=mvProbit(formula,data=dataset,start=startv,iterlim=100)
  
  #1035 Extract the estimates and standard errors
  estimates <- coef(result)
  standard_errors <- sqrt(diag(vcov(result)))
  
  # Combine the results into a data frame
  results_df <- data.frame(
    Parameter = parameters,
    Estimate = estimates,
    Std_Error = standard_errors
  )
  
  # Save the results to a CSV file for each dataset
  write.xlsx(results_df, paste0("C:/Users/u0118563/Documents/ml/r_", i, ".xlsx"), row.names = FALSE)
  
  # Print message that the dataset is done 
  print(paste("Dataset", i, "Done"))
}


# Start timer
write.xlsx(as.data.frame(rep(0,5)), paste0("C:/Users/u0118563/Documents/ml/START.xlsx"), row.names = FALSE)

# Loop over datasets simulated_0 to simulated_99
for(i in 0:99){
  while(length(dir(tempdir()))>100){
    # Remove temporary files 
    tmp_dir <- tempdir()
    files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    file.remove(files)
  }
print(length(dir(tempdir())))
analyze_dataset(i)
}

