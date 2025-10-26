rm(list=ls())
library(readxl)
library(dplyr)
library(scales)
library(reshape2)
# Setting work directory
setwd("C:\\Users\\SakshamMehta_k\\OneDrive - CHEORS\\Desktop\\Internal_Advanced_R_modelling_project")
getwd()

# Importing excel file
advances_R_excel=
  read_xlsx("C:\\Users\\SakshamMehta_k\\OneDrive - CHEORS\\Desktop\\Internal_Advanced_R_modelling_project\\R_final_docs\\R_PSM_Model_Input_Sheet - 26 September, 2025_10 AM (EST).xlsx",sheet = "Parameters")

age_disutilities=
  read_xlsx("C:\\Users\\SakshamMehta_k\\OneDrive - CHEORS\\Desktop\\Internal_Advanced_R_modelling_project\\R_final_docs\\R_PSM_Model_Input_Sheet - 26 September, 2025_10 AM (EST).xlsx",sheet = "Age related disutilities")

general_mortalities=
  read_xlsx("C:\\Users\\SakshamMehta_k\\OneDrive - CHEORS\\Desktop\\Internal_Advanced_R_modelling_project\\R_final_docs\\R_PSM_Model_Input_Sheet - 26 September, 2025_10 AM (EST).xlsx",sheet = "General Mortalities")

kmdata =
  read_xlsx("C:\\Users\\SakshamMehta_k\\OneDrive - CHEORS\\Desktop\\Internal_Advanced_R_modelling_project\\R_final_docs\\R_PSM_Model_Input_Sheet - 26 September, 2025_10 AM (EST).xlsx",
                    sheet = "Model Parameters")

# Approach selection, cutoff point storing, Median Duration storing
approach_selection <- as.character(as.matrix(kmdata[8:11, 3]))
distribution_selection <- as.character(as.matrix(kmdata[8:11, 4]))
cutoff_point <- as.numeric(as.matrix(kmdata[8:11, 7]))
med_cutoff_point <- as.numeric(as.matrix(kmdata[8:11, 8]))

data_transformed <- as.numeric(as.matrix(kmdata[16:145, 1:8]))
km_data <- matrix(data_transformed, nrow = 130)


# Assigning values
run_model <- function(inputs_df) {
  for(i in 1:nrow(inputs_df)) {
    assign(inputs_df$r_label[i], inputs_df$Mean[i], envir = .GlobalEnv)
  }
 

# Parameters
cycles <- 0:(Years*12)
advances_R_excel$Mean<-as.numeric(advances_R_excel$Mean) 
Years

# Transforming KM data as per the model 
km_data_transformed <- matrix(cycles, nrow = length(cycles))

DetermineKMdata=function(n){
  row_num1 = match(TRUE, km_data[, 1] > n) - 1
  row_num2 = match(TRUE, km_data[, 3] > n) - 1
  row_num3 = match(TRUE, km_data[, 5] > n) - 1
  row_num4 = match(TRUE, km_data[, 7] > n) - 1
  DetermineKMdata1 = km_data[row_num1, 2]
  DetermineKMdata2 = km_data[row_num2, 2]
  DetermineKMdata3 = km_data[row_num3, 2]
  DetermineKMdata4 = km_data[row_num4, 2]
  DetermineKMdata = matrix(c(DetermineKMdata1, DetermineKMdata2, DetermineKMdata3, DetermineKMdata4), nrow = 1)
  return(DetermineKMdata)
}

  data2 <- t(apply(km_data_transformed,1,DetermineKMdata))
  km_data_transformed <- cbind(km_data_transformed, data2)
  

# Function to create the survival curve for extrapolation using one-piece fitting 
# only using following distributions
# Exponential
# Weibull
# Lognormal
# Gompertz
# Loglogistics

surv_cur_generator <- function(cycles, distn, par_1, par_2){
  if(distn == "Exponential"){
    curve <- 1 - pexp(cycles, rate = par_1)
  }else if(distn == "Weibull"){
    curve <- 1 - pweibull(cycles, shape = par_1, scale = par_2)
  }else if(distn == "LogNormal"){
    curve <- 1 - plnorm(cycles, meanlog = par_1, sdlog = par_2)
  }else if(distn == "Gompertz"){
    curve <- 1 - pgompertz(cycles, shape = par_1, rate = par_2)
  }else if(distn == "LogLogistics"){
    curve <- 1 - plogis(cycles, location = par_1, scale = par_2)
  }
  return(curve)
}

# Function to create the survival curve for extrapolation using one-piece fitting, 
# two piece fitting, median treatment duration

# Approach 1 <- One-piece fitting
# Inputs needed <- Cycles, Distribution, Parameters

# Approach 2 <- Median duration
# Inputs needed <- Cycles, Treatment Duration

# Approach 3 <- Two-piece fitting
# Inputs needed <- Cycles, Cutoff-point, KM data, Distribution, Parameters
surv_cur_generator_global <- function(approach, cycles, cutoff, km_data, distn,
                                      par_1, par_2, med_dur){
  curve=c()
  
  if(approach == "One-piece"){
    curve <- surv_cur_generator(cycles, distn, par_1, par_2)
    
  }else if(approach == "Median Duration"){
    curve <- surv_cur_generator(cycles, "Exponential", log(2)/med_dur,) 
    
  }else if(approach == "Two-piece"){
    curve[1:(cutoff+1)] <- km_data[1:(cutoff+1), 2]
    
    onepiece <- surv_cur_generator(cycles, distn, par_1, par_2)
    adj_fac <- onepiece[-1]/onepiece[-length(onepiece)]
    
    last_pt <- km_data[cutoff+1, 2]
    
    rec_mul <- function(start, factors, n){
      if(n==1){
        return(c(start))
      }else{
        prev_vec <- rec_mul(start, factors, n-1)
        next_val <- prev_vec[length(prev_vec)]*factors[n-1]
        return(c(prev_vec,next_val))
      }
    }
    adjusted_curve <- rec_mul(adj_fac[1] * last_pt, adj_fac[2:length(adj_fac)], length(adj_fac))
    curve[(cutoff+2):(length(cycles)+1)] <- adjusted_curve[1:(length(cycles)+1-cutoff-2+1)]
  }
  return(curve)
} 

# Creating PF sheet matrix for Guma
guma_pf_matrix <- matrix(0, length(cycles), 35)
guma_pf_matrix[, 1] <- cycles
guma_pf_matrix[, 2] <- cycles/12

# Discount rate
guma_pf_matrix[, 3] <- 1/(1+con_DR_QALYs)^guma_pf_matrix[,2]
guma_pf_matrix[, 4] <- 1/(1+con_DR_costs)^guma_pf_matrix[,2]

# Creating PF sheet matrix for Savo
savo_pf_matrix <- matrix(0, length(cycles), 35)
savo_pf_matrix[, 1] <- cycles
savo_pf_matrix[, 2] <- cycles/12

# Discount rate
savo_pf_matrix[, 3] <- 1/(1+con_DR_QALYs)^savo_pf_matrix[,2]
savo_pf_matrix[, 4] <- 1/(1+con_DR_costs)^savo_pf_matrix[,2]

# # Survival Probabilities
# savo_pf_matrix[, 5] <- surv_cur_generator_global(approach_selection[3], cycles, 
#                                                  cutoff_point[3], km_data_transformed[,4], 
#                                                  distribution_selection[3],
#                                                  savo_OS_exp_lambda, , )
# savo_pf_matrix[, 6] <- surv_cur_generator_global(approach_selection[4], cycles,
#                                                  cutoff_point[4], km_data_transformed[,5],
#                                                  distribution_selection[4],
#                                                  savo_PFS_lognorm_lambda, savo_PFS_lognorm_sigma, )

# # Survival probabilities
# apply_general_mortalities <- function(matrix,OS_prob=c(approach_select, cutoff_poi, km_data_tra, distribution_sel, guma_OS_exp_lambda,savo_OS_exp_lambda)){
#   if (general_mortality_code==1){
#     matrix[,5] <- pmin(surv_cur_generator_global(approach_select, cycles, cutoff_poi, km_data_tra, distribution_sel,
#                                                  guma_OS_exp_lambda, savo_OS_exp_lambda, ),general_mortalities$`Survival probability`[1:109])
#   } 
#   else {
#     matrix[,5] <- surv_cur_generator_global(approach_select, cycles, cutoff_poi, km_data_tra, distribution_sel,
#                                             guma_OS_exp_lambda, savo_OS_exp_lambda, )
#   }
# }
# 
# guma_pf_matrix[, 5] <- apply_general_mortalities(approach_selection[1], cycles, 
#                                                  cutoff_point[1], km_data_transformed[,2],
#                                                  distribution_selection[1], guma_OS_exp_lambda, )

if(general_mortality_code==1){
  guma_pf_matrix[, 5] <- pmin(surv_cur_generator_global(approach_selection[1], 
                                                              cycles, cutoff_point[1], km_data_transformed[,2], 
                                                              distribution_selection[1],
                                                              guma_OS_exp_lambda, , ), 
                                                              general_mortalities$`Survival probability`[1:109])
}else{
  guma_pf_matrix[, 5] <- surv_cur_generator_global(approach_selection[1], 
                                                              cycles, cutoff_point[1], km_data_transformed[,2], 
                                                              distribution_selection[1],
                                                              guma_OS_exp_lambda,, )
}
guma_pf_matrix[, 6] <- surv_cur_generator_global(approach_selection[2], cycles, cutoff_point[2], 
                                                 km_data_transformed[,3], distribution_selection[2],
                                                 guma_PFS_lognorm_lambda, guma_PFS_lognorm_sigma, )
#Treatment waning
apply_treatment_waning <- function(guma, savo, col_OS = 5, col_PFS = 6, start_year, end_year, apply = TRUE) {
  if (!apply) return(guma)
  
  s <- start_year * 12; e <- end_year * 12
  r_OS <- ifelse(is.finite(guma[e+1, col_OS]) && guma[e+1, col_OS] != 0, savo[e+1, col_OS] / guma[e+1, col_OS], 1)
  r_PFS <- ifelse(is.finite(guma[e+1, col_PFS]) && guma[e+1, col_PFS] != 0, savo[e+1, col_PFS] / guma[e+1, col_PFS], 1)
  
  f_OS <- c(rep(1, s), seq(1, r_OS, length.out = e - s + 1), rep(r_OS, nrow(guma) - e - 1))
  f_PFS <- c(rep(1, s), seq(1, r_PFS, length.out = e - s + 1), rep(r_PFS, nrow(guma) - e - 1))
  
  guma[, col_OS] <- guma[, col_OS] * f_OS
  guma[, col_PFS] <- guma[, col_PFS] * f_PFS
  return(guma)
}
#Applying treatment waning based on selection
guma_pf_matrix <- apply_treatment_waning(guma = guma_pf_matrix, savo = savo_pf_matrix,
                                         start_year = TW_start_year,
                                         end_year = TW_end_year,
                                         apply = (con_TW == 1)
)
guma_pf_matrix[, 7] <- pmin(guma_pf_matrix[, 5], guma_pf_matrix[, 6])
guma_pf_matrix[, 8] <- pmax(guma_pf_matrix[, 5]-guma_pf_matrix[, 6], 0)
guma_pf_matrix[, 9] <- 1 - guma_pf_matrix[, 5]
guma_pf_matrix[, 10] <- 
  (guma_pf_matrix[, 7] + guma_pf_matrix[, 8] + guma_pf_matrix[, 9])==1

#Half cycle correction
apply_hcc <- function(matrix, source_col, target_col) {
  if (con_HCC == 1) {
    matrix[1, target_col] <- matrix[1, source_col]
    matrix[2:nrow(matrix), target_col] <- (matrix[1:(nrow(matrix)-1), source_col] + matrix[2:nrow(matrix), source_col]) / 2
  } else {
    matrix[, target_col] <- matrix[, source_col]
  }
  
  return(matrix)
}

#Applying HCC in PF health state
guma_pf_matrix <- apply_hcc(guma_pf_matrix, 7, 11)
guma_pf_matrix[, 11]

#Applying HCC in PD health state
guma_pf_matrix <- apply_hcc(guma_pf_matrix, 8, 12)
guma_pf_matrix[, 12]

#Applying HCC in Death health state
guma_pf_matrix <- apply_hcc(guma_pf_matrix, 9, 13)
guma_pf_matrix[, 13]

guma_pf_matrix[, 14] <- 
  (guma_pf_matrix[, 11] + guma_pf_matrix[, 12] + guma_pf_matrix[, 13])==1


# Calculating LYs
guma_pf_matrix[, 15] <- guma_pf_matrix[, 11]/12
guma_pf_matrix[, 16] <- guma_pf_matrix[, 12]/12


# Calculating QALYs
guma_pf_matrix[, 17] <- guma_pf_matrix[, 15]*util_PFS
guma_pf_matrix[, 18] <- guma_pf_matrix[, 16]*util_PD
# One off AE disutilities
apply_ae_disutility <- function(matrix,x,y){
    if (ae_disutility_code==1){
    matrix[,19] <- c(sum(advances_R_excel[50:60,2]*advances_R_excel[x:y,2]),rep(0,108))
  } 
  else {
    matrix[,19] <- (sum(advances_R_excel[50:60,2]*advances_R_excel[x:y,2])*matrix[,11])/nrow(matrix)
  }
}
guma_pf_matrix[,19] <- apply_ae_disutility(guma_pf_matrix,61,71)

#Age related disutilities
apply_age_disutility <- function(matrix){
  if (age_related_disutility==1){
    matrix[,20] <- ((1-matrix[,13])*age_disutilities$`Age-related disutility`[1:109])
  } 
  else {
    matrix[,20] <- 0
  }
}
guma_pf_matrix[,20] <- apply_age_disutility(guma_pf_matrix)

# Undiscounted costs per cycle
# Drug cost
guma_pf_matrix[, 21] <- guma_per_cycle_drug_cost*guma_pf_matrix[,11]
guma_pf_matrix[, 21]

# One off AE cost
apply_ae_cost <- function(matrix,x,y){
  if (ae_cost_code==1){
    matrix[,22] <- c(sum(advances_R_excel[x:y,2]*advances_R_excel[39:49,2]),rep(0,108))
  } 
  else {
    matrix[,22] <- (sum(advances_R_excel[x:y,2]*advances_R_excel[39:49,2])*matrix[,11])/nrow(guma_pf_matrix)
  }
}
guma_pf_matrix[,22] <- apply_ae_cost(guma_pf_matrix,61,71)
sum(guma_pf_matrix[, 22])

# Subsequent treatment cost
guma_pf_matrix[, 23] <- guma_pf_matrix[, 12]*peme_cisp_per_cycle_drug_cost
guma_pf_matrix[, 23]
# Disease Management first cycle cost
guma_pf_matrix[1, 24] <- sum(advances_R_excel[18:26,2])
# Disease Management subsequent cycles cost
guma_pf_matrix[2:109,24] <- ((sum(advances_R_excel[27:31,2])/42*(365.25/12))*
                               guma_pf_matrix[,11][-1]) +
  ((sum(advances_R_excel[32:38,2])/21*(365.25/12))*guma_pf_matrix[,12][-1])
guma_pf_matrix[,24]

# Terminal care cost
guma_pf_matrix[, 25] <- (guma_pf_matrix[, 13]-lag(guma_pf_matrix[,13],default = 0))*
  terminal_cost
guma_pf_matrix[,25]

# Discounted QALYs
guma_pf_matrix[, 26:30] <- guma_pf_matrix[, 15:19]*guma_pf_matrix[, 3]

#Discounted Costs
guma_pf_matrix[, 31:35] <- guma_pf_matrix[, 21:25]*guma_pf_matrix[, 4]


# Assigning column names
colnames(guma_pf_matrix)[1:35]=
  c("Cycles","Years","discounting factor of utilities","discounting factor of costs",
    "survival probabilities of OS",
    "survival probabilities of PFS","health state occupancy of PFS",
    "health state occupancy of PD","health state occupancy of Death","Check",
    "Half-cycle correction (PFS)","Half-cycle correction (PD)","Half-cycle correction (Death)",
    "Check","undiscounted LYs of PFS","undiscounted LYs of PD","undiscounted QALYs of PFS",
    "undiscounted QALYs of PD","AE disutility","Age-related disutility","undiscounted drug cost ($)",
    "AE cost ($)","undiscounted subsequent cost ($)",
    "undiscounted disease management cost ($)","undiscounted terminal cost ($)",
    "discounted LYs of PFS","discounted LYs of PD","discounted QALYs of PFS",
    "discounted QALYs of PD","AE disutility","discounted drug cost ($)",
    "AE cost ($)","discounted subsequent cost ($)",
    "discounted disease management cost ($)","discounted terminal cost ($)")



# Survival Probabilities
if(general_mortality_code==1){
  savo_pf_matrix[, 5] <- pmin(surv_cur_generator_global(approach_selection[3],
                                              cycles, cutoff_point[3], km_data_transformed[,6], 
                                              distribution_selection[3],
                                              savo_OS_exp_lambda, , ), 
                                              general_mortalities$`Survival probability`[1:109])
  }else{
    savo_pf_matrix[, 5] <- surv_cur_generator_global(approach_selection[3],
                                              cycles, cutoff_point[3], km_data_transformed[,6], 
                                              distribution_selection[3],
                                              savo_OS_exp_lambda,, )
  }
savo_pf_matrix[, 6] <- surv_cur_generator_global(approach_selection[4], cycles, cutoff_point[4], km_data_transformed[,8], distribution_selection[4],
                                                 savo_PFS_lognorm_lambda, savo_PFS_lognorm_sigma, )


savo_pf_matrix[, 7] <- pmin(savo_pf_matrix[, 5], savo_pf_matrix[, 6])
savo_pf_matrix[, 8] <- pmax(savo_pf_matrix[, 5]-savo_pf_matrix[, 6], 0)
savo_pf_matrix[, 9] <- 1 - savo_pf_matrix[, 5]
savo_pf_matrix[, 10] <- 
  (savo_pf_matrix[, 7] + savo_pf_matrix[, 8] + savo_pf_matrix[, 9])==1

#Applying HCC in PF health state
savo_pf_matrix <- apply_hcc(savo_pf_matrix, 7, 11)
savo_pf_matrix[, 11]

#Applying HCC in PD health state
savo_pf_matrix <- apply_hcc(savo_pf_matrix, 8, 12)
savo_pf_matrix[, 12]

#Applying HCC in Death health state
savo_pf_matrix <- apply_hcc(savo_pf_matrix, 9, 13)
savo_pf_matrix[, 13]

savo_pf_matrix[, 14] <- 
  (savo_pf_matrix[, 11] + savo_pf_matrix[, 12] + savo_pf_matrix[, 13])==1

savo_pf_matrix[, 14] <- 
  (savo_pf_matrix[, 11] + savo_pf_matrix[, 12] + savo_pf_matrix[, 13])==1


# Calculating LYs
savo_pf_matrix[, 15] <- savo_pf_matrix[, 11]/12
savo_pf_matrix[, 16] <- savo_pf_matrix[, 12]/12

# Calculating QALYs
savo_pf_matrix[, 17] <- savo_pf_matrix[, 15]*util_PFS
savo_pf_matrix[, 18] <- savo_pf_matrix[, 16]*util_PD
savo_pf_matrix[,19] <- apply_ae_disutility(savo_pf_matrix,72,82)
savo_pf_matrix[,20] <- apply_age_disutility(savo_pf_matrix)

#undiscounted drug cost
# Drug cost
savo_pf_matrix[, 21] <- savo_per_cycle_drug_cost*savo_pf_matrix[,11]
# One off AE cost
savo_pf_matrix[,22] <- apply_ae_cost(savo_pf_matrix,72,82)
# Subsequent treatment cost
savo_pf_matrix[, 23] <- savo_pf_matrix[, 12]*peme_cisp_per_cycle_drug_cost
# Disease management first cycle cost
savo_pf_matrix[1, 24] <- sum(advances_R_excel[18:26,2])
# Disease management subsequent cycles cost
savo_pf_matrix[2:109,24] <- ((sum(advances_R_excel[27:31,2])/42*(365.25/12))*
                               savo_pf_matrix[,11][-1]) +
  ((sum(advances_R_excel[32:38,2])/21*(365.25/12))*savo_pf_matrix[,12][-1])
# Terminal care cost
savo_pf_matrix[, 25] <- (savo_pf_matrix[, 13]-lag(savo_pf_matrix[,13],default = 0))*
  terminal_cost
savo_pf_matrix[,25]

# Discounted costs
savo_pf_matrix[, 26:30] <- savo_pf_matrix[, 15:19]*savo_pf_matrix[, 3]
savo_pf_matrix[, 31:35] <- savo_pf_matrix[, 21:25]*savo_pf_matrix[, 4]


# Assigning column names
colnames(savo_pf_matrix)[1:35]=
  c("Cycles","Years","discounting factor of utilities","discounting factor of costs",
    "survival probabilities of OS",
    "survival probabilities of PFS","health state occupancy of PFS",
    "health state occupancy of PD","health state occupancy of Death","Check",
    "Half-cycle correction (PFS)","Half-cycle correction (PD)","Half-cycle correction (Death)",
    "Check","undiscounted LYs of PFS","undiscounted LYs of PD","undiscounted QALYs of PFS",
    "undiscounted QALYs of PD","AE disutility","Age-related disutility","undiscounted drug cost ($)",
    "AE cost ($)","undiscounted subsequent cost ($)",
    "undiscounted disease management cost ($)","undiscounted terminal cost ($)",
    "discounted LYs of PFS","discounted LYs of PD","discounted QALYs of PFS",
    "discounted QALYs of PD","AE disutility","discounted drug cost ($)",
    "AE cost ($)","discounted subsequent cost ($)",
    "discounted disease management cost ($)","discounted terminal cost ($)")



# Costs
# Table for disaggregated total costs
guma_final_df=as.data.frame(guma_pf_matrix)


savo_final_df=as.data.frame(savo_pf_matrix)


disaggregated_dis_costs=data.frame(
  "Cost_category"=c("Total drug cost","AE Management Cost","Subsequent treatment cost",
                    "Disease Management Cost","Terminal care cost"),
  "Gumarontinib ($)"=comma(c(sum(guma_final_df$`discounted drug cost ($)`),
                             sum(guma_final_df$`AE cost ($)`),
                             sum(guma_final_df$`discounted subsequent cost ($)`),
                             sum(guma_final_df$`discounted disease management cost ($)`),
                             sum(guma_final_df$`discounted terminal cost ($)`)),
                           accuracy = .01),
  "Savolitinib ($)"=comma(c(sum(savo_final_df$`discounted drug cost ($)`),
                            sum(savo_final_df$`AE cost ($)`),
                            sum(savo_final_df$`discounted subsequent cost ($)`),
                            sum(savo_final_df$`discounted disease management cost ($)`),
                            sum(savo_final_df$`discounted terminal cost ($)`)),
                          accuracy = 0.01)
)

final_disaggregated_dis_costs=
  rbind(disaggregated_dis_costs,
        c("Total discounted costs",
          comma(sum(sum(guma_final_df$`discounted drug cost ($)`),
                    sum(guma_final_df$`AE cost ($)`),
                    sum(guma_final_df$`discounted subsequent cost ($)`),
                    sum(guma_final_df$`discounted disease management cost ($)`),
                    sum(guma_final_df$`discounted terminal cost ($)`)),
                accuracy = .01),
          comma(sum(sum(savo_final_df$`discounted drug cost ($)`),
                    sum(savo_final_df$`AE cost ($)`),
                    sum(savo_final_df$`discounted subsequent cost ($)`),
                    sum(savo_final_df$`discounted disease management cost ($)`),
                    sum(savo_final_df$`discounted terminal cost ($)`)),
                accuracy = .01)))

final_disaggregated_dis_costs$Gumarontinib....=
  format(final_disaggregated_dis_costs$Gumarontinib....,big.mark = ",")
final_disaggregated_dis_costs$Gumarontinib....



# Life Years
# Table of disaggregated life years
disaggregated_dis_lys=data.frame(
  "Health_State"=c("PF","PD","Total"),
  "Gumarontinib"=
    comma(c(sum(guma_final_df$`discounted LYs of PFS`),
            sum(guma_final_df$`discounted LYs of PD`),
            sum(sum(guma_final_df$`discounted LYs of PFS`),
                sum(guma_final_df$`discounted LYs of PD`))),accuracy = .001),
  "Savolitinib"=
    comma(c(sum(savo_final_df$`discounted LYs of PFS`),
            sum(savo_final_df$`discounted LYs of PD`),
            sum(sum(savo_final_df$`discounted LYs of PFS`),
                sum(savo_final_df$`discounted LYs of PD`))),accuracy = .001))




# Quality Adjusted Life Years
# Table of disaggregated QALYs
disaggregated_dis_qalys=data.frame(
  "Health_State"=c("PF","PD","AE Disutility","Age-related disutility","Total"),
  "Gumarontinib"=
    comma(c(sum(guma_final_df$`discounted QALYs of PFS`),
            sum(guma_final_df$`discounted QALYs of PD`),
            sum(guma_final_df$`AE disutility`),
            sum(guma_final_df$`Age-related disutility`),
            sum(sum(guma_final_df$`discounted QALYs of PFS`),
                sum(guma_final_df$`discounted QALYs of PD`),
                sum(guma_final_df$`Age-related disutility`),
                sum(guma_final_df$`AE disutility`))),accuracy = .001),
  "Savolitinib"=
    comma(c(sum(savo_final_df$`discounted QALYs of PFS`),
            sum(savo_final_df$`discounted QALYs of PD`),
            sum(savo_final_df$`AE disutility`),
            sum(savo_final_df$`Age-related disutility`),
            sum(sum(savo_final_df$`discounted QALYs of PFS`),
                sum(savo_final_df$`discounted QALYs of PD`),
                sum(savo_final_df$`Age-related disutility`),
                sum(savo_final_df$`AE disutility`))),accuracy = .001))




# Calculating all the total life years, qalys and costs
guma_total_dis_lys=
  sum(sum(guma_final_df$`discounted LYs of PFS`),
      sum(guma_final_df$`discounted LYs of PD`))  
savo_total_dis_lys=
  sum(sum(savo_final_df$`discounted LYs of PFS`),
      sum(savo_final_df$`discounted LYs of PD`))
guma_total_dis_qalys=
  sum(sum(guma_final_df$`discounted QALYs of PFS`),
      sum(guma_final_df$`discounted QALYs of PD`),
      sum(guma_final_df$`Age-related disutility`),
      sum(guma_final_df$`AE disutility`))
savo_total_dis_qalys=
  sum(sum(savo_final_df$`discounted QALYs of PFS`),
      sum(savo_final_df$`discounted QALYs of PD`),
      sum(savo_final_df$`Age-related disutility`),
      sum(savo_final_df$`AE disutility`))
guma_total_dis_costs=
  sum(sum(guma_final_df$`discounted drug cost ($)`),
      sum(guma_final_df$`AE cost ($)`),
      sum(guma_final_df$`discounted subsequent cost ($)`),
      sum(guma_final_df$`discounted disease management cost ($)`),
      sum(guma_final_df$`discounted terminal cost ($)`))
savo_total_dis_costs=
  sum(sum(savo_final_df$`discounted drug cost ($)`),
      sum(savo_final_df$`AE cost ($)`),
      sum(savo_final_df$`discounted subsequent cost ($)`),
      sum(savo_final_df$`discounted disease management cost ($)`),
      sum(savo_final_df$`discounted terminal cost ($)`))
incr_costs=(guma_total_dis_costs-savo_total_dis_costs)
incr_lys=(guma_total_dis_lys-savo_total_dis_lys)
incr_qalys=(guma_total_dis_qalys-savo_total_dis_qalys)
icer_by_ly=incr_costs/incr_lys
icer_by_qaly=incr_costs/incr_qalys


# Final ICER Table
final_df=data.frame(
  "Group"=c("Gumarontinib","Savolitinib"),
  "Discounted_LYs"=comma(c(guma_total_dis_lys,savo_total_dis_lys),accuracy=0.01),
  "Discounted_QALYs"=comma(c(guma_total_dis_qalys,savo_total_dis_qalys),accuracy=0.01),
  "Total_discounted_cost_($)"=comma(c(guma_total_dis_costs,savo_total_dis_costs),accuracy=0.01),
  "Incremental_LYs"=c("",comma(incr_lys,accuracy=0.01)),
  "Incremental_QALYs"=c("",comma(incr_qalys,accuracy=0.01)),
  "Incremental_cost_($)"=c("",comma(incr_costs,accuracy=0.01)),
  "ICER_by_LY"=c("",comma(icer_by_ly,accuracy=0.01)),
  "ICER_by_QALY"=c("",comma(icer_by_qaly,accuracy=0.01))
)

NMB = incr_qalys*con_WTP - incr_costs
# Step 3: Return the results
return(list(guma_pf_matrix = guma_pf_matrix,
            savo_pf_matrix = savo_pf_matrix,
            guma_final_df = guma_final_df,
            savo_final_df = savo_final_df,
            final_disaggregated_dis_costs = final_disaggregated_dis_costs,
            disaggregated_dis_lys = disaggregated_dis_lys,
            disaggregated_dis_qalys = disaggregated_dis_qalys,
            final_df = final_df,
            guma_total_dis_costs = guma_total_dis_costs,
            savo_total_dis_costs = savo_total_dis_costs,
            guma_total_dis_qalys = guma_total_dis_qalys,
            savo_total_dis_qalys = savo_total_dis_qalys,
            incr_costs = incr_costs, incr_qalys = incr_qalys, 
            icer_by_qaly = icer_by_qaly,
            NMB = NMB))
}
base_case_results<- run_model(advances_R_excel)
View(base_case_results$guma_final_df)
View(base_case_results$savo_final_df)
View(base_case_results$final_disaggregated_dis_costs)
View(base_case_results$disaggregated_dis_lys)
View(base_case_results$disaggregated_dis_qalys)
View(base_case_results$final_df)

#setting OWSA ICER
owsa_results <- data.frame(Parameter = character(), Low = numeric(), High = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(advances_R_excel)) {
  param_name <- advances_R_excel$r_label[i]
  lower_val <- advances_R_excel$Lower[i]
  upper_val <- advances_R_excel$Upper[i]
  # Skip if either bound is missing (NA or blank)
  if (is.na(lower_val) || is.na(upper_val) || lower_val == "" || upper_val == "") {
    next
  }
  temp_input <- advances_R_excel
  # Lower bound
  temp_input$Mean[i] <- lower_val
  temp_input$Mean <- as.numeric(temp_input$Mean)  
  low_icer <- run_model(temp_input)$icer_by_qaly
  # Upper bound
  temp_input$Mean[i] <- upper_val
  temp_input$Mean <- as.numeric(temp_input$Mean)
  high_icer <- run_model(temp_input)$icer_by_qaly
  owsa_results <- rbind(owsa_results, data.frame(
    Parameter = param_name,
    Low = low_icer,
    High = high_icer
  ))
}
View(owsa_results)
library(ggplot2)
#format ICER
icer_val<- paste0("$", format(round(base_case_results$icer_by_qaly, 0), big.mark = ","))

#Calculate Impact
owsa_results$Impact<- pmax(abs(owsa_results$Low - base_case_results$icer_by_qaly),
                           abs(owsa_results$High - base_case_results$icer_by_qaly))
#Top 15 parameters
owsa_top15<- owsa_results[order(owsa_results$Impact, decreasing = TRUE), ][1:15, ]
#Reordering
owsa_top15$Parameter<- factor(owsa_top15$Parameter, 
                              levels = owsa_top15$Parameter[order(owsa_top15$Impact, decreasing = TRUE)])

#Plot Tornado
ggplot() +
  geom_segment(data = owsa_top15, aes(x = Low, xend= base_case_results$icer_by_qaly, y= Parameter, yend = Parameter,
                                      color ="Lower Bound"), size = 4)+
  geom_segment(data = owsa_top15, aes(x = base_case_results$icer_by_qaly, xend = High, y= Parameter, yend = Parameter,
                                      color = "Upper Bound"), size = 4) +
  geom_vline(xintercept = base_case_results$icer_by_qaly, color = "black", size = 0.5) +
  scale_color_manual(name = "Bound", values = c("Lower Bound" = "red", "Upper Bound" = "blue"))+
  scale_y_discrete(limits = rev(levels(owsa_top15$Parameter)))+
  labs(x= paste0("ICER = ", icer_val), y= "", title = "Tornado Diagram (ICER)") +
  theme_minimal() +
  theme(legend.position = "right")


#setting OWSA NMB
owsa_results_NMB <- data.frame(Parameter = character(), Low = numeric(), High = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(advances_R_excel)) {
  param_name <- advances_R_excel$r_label[i]
  lower_val <- advances_R_excel$Lower[i]
  upper_val <- advances_R_excel$Upper[i]
  # Skip if either bound is missing (NA or blank)
  if (is.na(lower_val) || is.na(upper_val) || lower_val == "" || upper_val == "") {
    next
  }
  temp_input <- advances_R_excel
  # Lower bound
  temp_input$Mean[i] <- lower_val
  temp_input$Mean <- as.numeric(temp_input$Mean)  
  low_NMB <- run_model(temp_input)$NMB
  # Upper bound
  temp_input$Mean[i] <- upper_val
  temp_input$Mean <- as.numeric(temp_input$Mean)
  high_NMB <- run_model(temp_input)$NMB
  owsa_results_NMB <- rbind(owsa_results_NMB, data.frame(
    Parameter = param_name,
    Low = low_NMB,
    High = high_NMB
  ))
}
View(owsa_results_NMB)

#format NMB
NMB_val<- paste0("$", format(round(base_case_results$NMB, 0), big.mark = ","))

#Calculate Impact
owsa_results_NMB$Impact<- pmax(abs(owsa_results_NMB$Low - base_case_results$NMB),
                               abs(owsa_results_NMB$High - base_case_results$NMB))
#Top 15 parameters
owsa_top15_NMB<- owsa_results_NMB[order(owsa_results_NMB$Impact, 
                                        decreasing = TRUE), ][1:15, ]
#Reordering
owsa_top15_NMB$Parameter<- factor(owsa_top15_NMB$Parameter, 
                                  levels = owsa_top15_NMB$Parameter[order(owsa_top15_NMB$Impact, decreasing = TRUE)])

#Plot Tornado NMB
ggplot() +
  geom_segment(data = owsa_top15_NMB, aes(x = Low, xend= base_case_results$NMB, y= Parameter, yend = Parameter,
                                          color ="Lower Bound"), size = 4)+
  geom_segment(data = owsa_top15_NMB, aes(x = base_case_results$NMB, xend = High, y= Parameter, yend = Parameter,
                                          color = "Upper Bound"), size = 4) +
  geom_vline(xintercept = base_case_results$NMB, color = "black", size = 0.5) +
  scale_color_manual(name = "Bound", values = c("Lower Bound" = "red", "Upper Bound" = "blue"))+
  scale_y_discrete(limits = rev(levels(owsa_top15_NMB$Parameter)))+
  labs(x= paste0("NMB = ", NMB_val), y= "", title = "Tornado Diagram (NMB)") +
  theme_minimal() +
  theme(legend.position = "right")

set.seed(1234)
n_iter <- 1500

start_time <- Sys.time()  # Start time

# Function for generating random rows
sample_param <- function(mean, dist_type, n){
  if(is.na(dist_type) || dist_type == ""){
    return(rep(mean, n))
  }
  if(dist_type == "Gamma"){
    sd <- 0.1 * mean
    shape <- (mean / sd)^2
    scale <- sd^2 / mean
    return(rgamma(n, shape = shape, scale = scale))
  }
  if(dist_type == "Beta"){
    sd <- 0.1 * mean
    alpha <- ((1 - mean) / sd^2 - 1 / mean) * mean^2
    beta <- alpha * (1 / mean - 1)
    return(rbeta(n, alpha, beta))
  }
  stop("Unknown distribution type: ", dist_type)
}

# Generate PSA inputs iterations
psa_iterations <- data.frame(iter = 1:n_iter)
for(i in 1:nrow(advances_R_excel)){
  pname <- advances_R_excel$r_label[i]
  meanv <- advances_R_excel$Mean[i]
  distn <- advances_R_excel$Distribution[i]
  psa_iterations[[pname]] <- sample_param(meanv, distn, n_iter)
}

# PSA ICER calculation
param_cols <- setdiff(names(psa_iterations), "iter")
icer_vec <- apply(psa_iterations[, param_cols, drop = FALSE], 1, function(row_vals){
  temp_inp <- advances_R_excel
  for(j in seq_along(param_cols)){
    param_name <- param_cols[j]
    if(param_name %in% temp_inp$r_label){
      temp_inp$Mean[temp_inp$r_label == param_name] <- row_vals[j]
    }
  }
  run_model(temp_inp)
})

psa_results<- do.call(rbind, icer_vec)
psa_results<- as.data.frame(psa_results)
psa_results$guma_total_dis_costs<- as.numeric(psa_results$guma_total_dis_costs)
psa_results$savo_total_dis_costs<- as.numeric(psa_results$savo_total_dis_costs)
psa_results$guma_total_dis_qalys<- as.numeric(psa_results$guma_total_dis_qalys)
psa_results$savo_total_dis_qalys<- as.numeric(psa_results$savo_total_dis_qalys)
psa_results$incr_qalys<- as.numeric(psa_results$incr_qalys)
psa_results$incr_costs<- as.numeric(psa_results$incr_costs)
psa_results$icer_by_qaly<- as.numeric(psa_results$icer_by_qaly)

# Summarize
PSA_summary_table<- data.frame(
  Strategy = c("Comparator_(Savo)", "Intervention_(Guma)", "Incremental"),
  Mean_Costs = c(mean(psa_results$savo_total_dis_costs),
                 mean(psa_results$guma_total_dis_costs),
                 mean(psa_results$incr_costs)),
  Mean_QALYs = c(mean(psa_results$savo_total_dis_qalys),
                 mean(psa_results$guma_total_dis_qalys),
                 mean(psa_results$incr_qalys)),
  PSA_ICER = c(NA, NA, mean(psa_results$icer_by_qaly))
  
)
end_time <- Sys.time()  # End time
run_time <- end_time - start_time

# Output
View(PSA_summary_table)
print(paste("Run time:", run_time))

#CE Plane
base_incr_qaly<- base_case_results$incr_qalys
base_incr_cost<- base_case_results$incr_costs
mean_incr_qaly<- mean(psa_results$incr_qalys)
mean_incr_cost<- mean(psa_results$incr_costs)

# Make data frame for special ICER points
special_points <- data.frame(
  x = c(base_incr_qaly, mean_incr_qaly),
  y = c(base_incr_cost, mean_incr_cost),
  type = c("Base ICER", "Mean PSA ICER")
)

# Dummy data for WTP line (so it shows in legend)
wtp_line <- data.frame(type = "WTP Line")

ggplot() +
  # Iteration cloud
  geom_point(data = psa_results,
             aes(x = incr_qalys, y = incr_costs, color = "Iterations"), shape = 16,
             alpha = 0.4) +
  
  # Base + Mean ICER points
  geom_point(data = special_points,
             aes(x = x, y = y, color = type, shape = type),
             size = 4) +
  
  # WTP line (solid now)
  geom_abline(data = wtp_line,
              aes(slope = con_WTP, intercept = 0, color = "WTP Line", linetype = "WTP Line"),
              size = 1) +
  
  # Reference axes
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  labs(
    title = "Cost Effectiveness Plane",
    x = "Incremental QALYs",
    y = "Incremental Costs",
    
  ) +
  scale_color_manual(name = "Legend", 
                     limits =  c("Iterations", "Base ICER", "Mean PSA ICER",
                                 "WTP Line"),
                     values = c("Iterations" = "blue",
                                "Base ICER" = "green",
                                "Mean PSA ICER" = "red",
                                "WTP Line" = "darkorange")) +
  #Hiding extra legends
  scale_shape_manual(values = c("Base ICER" = 17, "Mean PSA ICER" = 18), 
                     guide = "none") +
  scale_linetype_manual(values = c("WTP Line" = "solid"),
                        guide = "none") +
  guides(
    color = guide_legend(override.aes = list(
      shape = c(16, 17, 18, NA),
      linetype = c(0,0,0,1),
      size = c(2.5, 4, 4,1)
    ))
  )+
  theme_minimal()+
  theme(legend.position = "right")


#CEAC
#WTP Thresholds
wtp_range<- seq(0, 150000, by = 10000)

#Store Prob
ceac_list<- lapply(wtp_range, function(wtp){
  #NMB
  nmb_inv<- psa_results$guma_total_dis_qalys*wtp - psa_results$guma_total_dis_costs
  nmb_comp<- psa_results$savo_total_dis_qalys*wtp - psa_results$savo_total_dis_costs
  
  #Check for which Strategy better
  prob_inv<- mean(nmb_inv > nmb_comp, na.rm = TRUE)
  prob_comp<- mean(nmb_comp > nmb_inv, na.rm = TRUE)
  
  data.frame(WTP = wtp, Guma = prob_inv, Savo = prob_comp)
})

#Combine results
ceac<- do.call(rbind, ceac_list)

#Reshape
ceac_long <- melt(ceac, id.vars = "WTP",
                  variable.name = "Strategy", value.name = "Probability")
#Plot CEAC
ggplot(ceac_long, aes(x = WTP, y = Probability, color = Strategy)) + 
  geom_line(size = 1.2) + 
  labs(
    title = "CEAC", x = "WTP_($_Per_QALY)", y = "Probability most CE"
  ) +
  theme_minimal() + 
  scale_color_manual(values = c("Guma" = "blue", "Savo" = "red"))
