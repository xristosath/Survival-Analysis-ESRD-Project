#load the necessary libraries
library(tidyverse)
library(survival)
set.seed(12)
#load data
esrd <- read.csv("esrd.csv", sep=";", header = T)

#half day in years
half_day_in_years=0.5/365.25

#check those with zero follow-up time
esrd %>%
  group_by(treatmnt)%>%
  summarise(
    n_zero_time = sum(dthtime == 0),
    n_total = n(),
    percent_zero = 100 * n_zero_time / n_total)

#prepare the data
esrd= esrd%>% mutate(
  dthtime=if_else(dthtime==0, half_day_in_years, dthtime) #we change those with death time zero to half day
)%>%
  mutate(
  #we create the survival object with death=1 as event indicator
  surv_obj= Surv(time=dthtime,event = dth),
  #factor all the categorical variables
  treatmnt= factor(treatmnt, levels = c(0,1), labels = c("Conservative","Haemodialysis")),
  sex= factor(sex, levels = c("F","M"), labels = c("Female", "Male")),
  dth=factor(dth,levels = c(0,1),labels = c("Censored","Died")),
  #we use across to factor all comorbidity columns
  across(copd:vascular, ~factor(.,levels = c(0,1), labels = c("No","Yes")))
  
  
)
# VERIFICATION CHECK
cat("Number of patients with time=0:", sum(esrd$dthtime == 0), "\n")
print(head(esrd$surv_obj))
##### 1st Question  ####
#Statistics

continuous_vars=c("age","dthtime") 
all_continuous_summaries= list()
for (variable_name in continuous_vars) {
  var_symbol= sym(variable_name) #covert the string "age" into variable that dplyr understands
  #get summary by tratment
  summary_by_group= esrd %>%
    group_by(treatmnt)%>%
    summarise(
      Mean = mean(!!var_symbol, na.rm = TRUE),
      SD = sd(!!var_symbol, na.rm = TRUE),
      Median = median(!!var_symbol, na.rm = TRUE),
      Q1 = quantile(!!var_symbol, 0.25, na.rm = TRUE),
      Q3 = quantile(!!var_symbol, 0.75, na.rm = TRUE)
    )
  #summary total
  summary_total= esrd %>%
    summarise(
      Mean = mean(!!var_symbol, na.rm = TRUE),
      SD = sd(!!var_symbol, na.rm = TRUE),
      Median = median(!!var_symbol, na.rm = TRUE),
      Q1 = quantile(!!var_symbol, 0.25, na.rm = TRUE),
      Q3 = quantile(!!var_symbol, 0.75, na.rm = TRUE)
    ) %>%
      mutate(treatmnt="Total")
  #stack and store them
  all_continuous_summaries[[variable_name]]=bind_rows(summary_by_group,summary_total) %>%
    mutate(Variable= variable_name)%>%  #add a column to know which variable this is
    select(Variable, treatmnt, everything()) #reprder columns
}
continuous_table= bind_rows(all_continuous_summaries)
all_continuous_summaries

#For categorical variabels now
categorical_vars= c("sex","copd","dm","hypert","heart","liver","neoplasia","vascular","dth")
all_categorical_summaries=list()
for (variable_name in categorical_vars) {
  var_symbol=sym(variable_name)
  
  #count and percentage by treatment
  summary_by_group= esrd %>%
    group_by(treatmnt) %>%
    count(Level= !!var_symbol) %>% #we rename the variable to Level so we can stack them all
    mutate(Percent= n/sum(n)*100)
  
  #count and percents for total
  summary_total= esrd %>%
    count(Level= !!var_symbol) %>% #we rename the variable to Level so we can stack them all
    mutate(Percent= n/sum(n)*100) %>%
    mutate(treatmnt="Total")
  
  #stack the group and total summaries
  all_categorical_summaries[[variable_name]]=bind_rows(summary_by_group,summary_total) %>%
    mutate(Variable= variable_name) %>% #add column to track the variable
    mutate(Variable, treatmnt, Level, n , Percent) #reorder
    
  
}
categorical_table= bind_rows(all_categorical_summaries)
all_categorical_summaries

#check again those with zero follow-up time 
esrd %>%
  group_by(treatmnt)%>%
  summarise(
    n_zero_time = sum(dthtime == 0),
    n_total = n(),
    percent_zero = 100 * n_zero_time / n_total
  )
esrd %>% filter(treatmnt == "Conservative", dthtime == 0) %>% head()
esrd%>%filter(treatmnt=="Conservative",id==16)%>% head()
 
#### Question 2 OVERALL DEATH RISK####

crude_by_trt= esrd%>%
  group_by(treatmnt)%>%
  summarise(
    n=n(), #number of rows in each group
    deaths=sum(dth=="Died"),
    censored=sum(dth=="Censored"),
    cuminc=deaths/n
  )%>%
  mutate(cuminc_pct=100* cuminc)
crude_by_trt
#add the total row in the table
total_row=esrd%>%
  summarise(
    n=n(),
    deaths=sum(dth=="Died"),
    censored=sum(dth=="Censored"),
    cuminc=deaths/n
  )%>%
  mutate(treatmnt="Total",cuminc_pct=100*cuminc)
  
total_row
final_risk_table=bind_rows(crude_by_trt,total_row)
final_risk_table

#Significance test

contingency_table=matrix(
  c(crude_by_trt$deaths[2],crude_by_trt$censored[2],
    crude_by_trt$deaths[1],crude_by_trt$censored[1]),
  nrow = 2, byrow = T,
  dimnames = list(Treatment=c("Haemodialysis","Conservative"),
                  Outcome=c("Died","Censored"))
)
chisq_test=chisq.test(contingency_table,correct = F) #no need for yates correction
contingency_table
chisq_test
#Effectiveness Calcs
risk_cons=final_risk_table$cuminc[final_risk_table$treatmnt=="Conservative"]
risk_haem=final_risk_table$cuminc[final_risk_table$treatmnt=="Haemodialysis"]

RR=risk_haem/risk_cons
RRR=1-RR
RR
RRR

####Question 3:Overall time-adjusted mortality rates per person year####
#calculate rates and CI
rate_summary= esrd%>%
  group_by(treatmnt)%>%
  summarise(
    Total_Deaths=sum(dth=="Died"),
    Total_PY=sum(dthtime)
  )%>%
  rowwise() %>% #to apply the test row by row
  mutate(
    #calculate exact poisson CI for the rate
    test_res=list(poisson.test(Total_Deaths,T=Total_PY)),
    Rate_per_PY=Total_Deaths/Total_PY,
    Rate_per_100_PY=Rate_per_PY*100,
    CI_Lower=test_res$conf.int[1]*100,
    CI_Upper=test_res$conf.int[2]*100
  ) %>%
  select(-test_res) #remove the temporary list column
rate_summary

rate_table_df <- data.frame(
  Treatment = rate_summary$treatmnt,
  `Total Deaths` = rate_summary$Total_Deaths,
  `Total Person-Years` = round(rate_summary$Total_PY, 1),
  `Rate (per 100 PY)` = sprintf("%.1f", rate_summary$Rate_per_100_PY),
  `95% CI` = sprintf("[%.1f, %.1f]", rate_summary$CI_Lower, rate_summary$CI_Upper)
)
rate_table_df
#Significance Test (Rate Ratio)
#we have rates so we will proceed with a poisson test /we could have answered the Q3 with Poisson Regression
deaths_cons= rate_summary$Total_Deaths[1]
py_cons= rate_summary$Total_PY[1]  # The real person-years
deaths_haem= rate_summary$Total_Deaths[2]
py_haem= rate_summary$Total_PY[2] # The real person-years 

poisson_test <- poisson.test(
  x = c(deaths_cons, deaths_haem),
  T = c(py_cons, py_haem) # This is now correct
)
poisson_test #the rate of 5.57 here is rate(conservative)/rate(haemodialysis) and we need it all the way around

#Effectiveness
rate_cons <- rate_summary$Rate_per_PY[rate_summary$treatmnt == "Conservative"]
rate_haem <- rate_summary$Rate_per_PY[rate_summary$treatmnt == "Haemodialysis"]

Rate_Ratio<- rate_haem/ rate_cons
RRR_q3 <- 1 - Rate_Ratio
Rate_Ratio
RRR_q3

####Question 4 Kaplan-Meier Curves####

#we start by fitting the KM survival model
km_fit=survfit(surv_obj~treatmnt, data = esrd)
summary(km_fit)
km_fit
#plot the curves
plot(km_fit,mark.time = F,main="Kaplan-Meier Survival Estimates by Treatment Modality",
     col = c("blue","red"),
     lty=1:2, xlab = "Time in Years", ylab = "Survival Probability")
legend("topright",lty = 1:2,col = c("blue","red"), legend = c("Conservative","Haemodialysis"),bty = "n")

#logrank test
logrank_test=survdiff(surv_obj~treatmnt,data = esrd)
logrank_test
#compare the median survival between the two groups
quantile(survfit(surv_obj~treatmnt,data = esrd),probs = 0.5,conf.int = F)

#peto peto
wilcoxon_test=survdiff(surv_obj~treatmnt,data = esrd,rho = 1)
wilcoxon_test

####Question 5 ####
#(a) estimate the survival probability at specific times
#we are going to use the km model from the previous question
summary_5a=summary(km_fit,times = c(0.5,1,3,5))
# We create a clean data frame from the summary object
table_5a_df <- data.frame(
  Time = summary_5a$time,
  Treatment = summary_5a$strata,
  Survival_Prob = summary_5a$surv
) %>%
# Pivot the table to the format requested
  pivot_wider(names_from = Treatment, values_from = Survival_Prob)
table_5a_df

#5(b) conditional survival probabilities at specific times

#we subset for patients who survived past 6 months (0.5 years)
km.fit_5b= survfit(Surv(dthtime-0.5,dth=="Died")~treatmnt,data = esrd,subset = (dthtime>0.5)) #the dthtime-0.5 makes for those their new t=0
additional_times=c(0.5,1,1.5,2)
summary_cond5b= summary(km.fit_5b,times = additional_times)

# Create a clean data frame 
table_5b <- data.frame(
  Time = summary_cond5b$time,
  Treatment = summary_cond5b$strata,
  Survival_Prob = summary_cond5b$surv
) %>%
  pivot_wider(names_from = Treatment, values_from = Survival_Prob) %>%
  # Add a column to clarify what is being estimated
  mutate(
    `Additional Time (t)` = additional_times,
    `Estimating` = paste0("S(t | T > 0.5) at t = ", c("0.5yr", "1yr", "1.5yr", "2yr"))
  ) %>%
  select(Estimating, `Additional Time (t)`, everything(), -Time) # Re-order
table_5b
#Question 6 : univariate cox model for age

cox.age=coxph(surv_obj~age,data = esrd)
uni.cox_summary=summary(cox.age)
uni.cox_summary
univ.cox_table= data.frame(
  Variable="Age (per year)",
  `Hazard Ration (HR)`=uni.cox_summary$conf.int[1],
  `Lower CI`=uni.cox_summary$conf.int[3],
  `Upper CI`=uni.cox_summary$conf.int[4],
  `p-value`=format.pval(uni.cox_summary$coefficients[5],eps = 0.001)
)
univ.cox_table
#we also need to check linearity
# We plot Martingale Residuals against Age.
# If the relationship is linear, the smooth line should be roughly straight/horizontal.
# If it curves U-shape or S-shape, "Age" might need a transformation (like age^2).


plot(esrd$age, residuals(cox.age, type = "martingale"),
     xlab = "Age", ylab = "Martingale Residuals",
     main = "Linearity Check for Age",
     pch = 20, col = "darkgray")
lines(lowess(esrd$age, residuals(cox.age, type = "martingale")), col = "red", lwd = 2)

#Question 7:Univariate Cox for treatment
cox.trt= coxph(surv_obj~treatmnt, data = esrd)
cox.trt_sum=summary(cox.trt)
cox.trt_sum
#table
cox_trt_table <- data.frame(
  Variable = "Haemodialysis (vs. Conservative)",
  Hazard_Ratio = cox.trt_sum$conf.int[1],        
  CI = paste0("[", round(cox.trt_sum$conf.int[3], 3), ", ", round(cox.trt_sum$conf.int[4], 3), "]"),            
  P_value = format.pval(cox.trt_sum$coefficients[5], eps = 0.001)
)
cox_trt_table 

#Question 8: PH assumption checking

#1st graphical approach with cloglog S(t)
#we use again the km_fit model we created for treatment in question 4
#if you haven't run the whole code the command is:
#km_fit=survfit(surv_obj~treatmnt, data = esrd)
#and we plot the cloglog
plot(km_fit,mark.time = F,fun = "cloglog", main="Log-Minus-Log Plot for Treatment",
     col = c("blue","red"),
     lty = 1:2,
     xlab = "Log(Time)",
     ylab = "ln[-ln(Survival Probabilities)]",
)
legend("topleft",
       legend = c("Conservative","Haemodialysis"),
       col = c("blue","red"),
       lty = 1:2,
       bty = "n"
       )
#2nd approach comparing the raw KM survival curves with those predicted by a Cox model
#we are using the cox model for treatment by question 7
#if needed, the command is : cox.trt= coxph(surv_obj~treatmnt, data = esrd)
#we now generate the predicted survival curves from the cox model
cox.pred=survfit(cox.trt,newdata = data.frame(treatmnt=c("Conservative","Haemodialysis")))

#now we first plot the fitted curves
plot(cox.pred, 
     col = c("blue", "red"),
     lty = 2, # Dashed line for Predicted
     lwd = 2,
     mark.time = FALSE,
     xlab = "Time in Years", 
     ylab = "Survival Probability", 
     main = "Observed (KM) vs. Predicted (Cox) Survival Curves \nBy Treatment"
)
#add the raw curves
lines(km_fit, 
      col = c("green", "orange"), 
      lty = 1, # Solid line for Observed
      lwd = 2, 
      mark.time = FALSE
)
legend("topright", 
       legend = c("Observed (KM): Cons.", "Observed (KM): Haem.", 
                  "Predicted (Cox): Cons.", "Predicted (Cox): Haem."),
       col = c("blue", "red", "green", "orange"),
       lty = c(1, 1, 2, 2), # 1=Solid, 2=Dashed
       lwd = 2,
       bty = "n",
       cex = 0.8
)
#3rd graphical approach through Scaled Schoenfeld Residfual
#Grambch & Therneau test
ph.test=cox.zph(cox.trt)
ph.test
#Extract the results matrix
ph.results <- ph.test$table

# Create a clean data frame for the table
ph.table <- data.frame(
  Variable = rownames(ph.results),
  `Chi-squared` = ph.results[, "chisq"],
  `df` = ph.results[, "df"],
  `P-value` = format.pval(ph.results[, "p"], eps = 0.001) # Format p < 0.001
)

# Clean up row names (optional, makes it look nicer)
ph.table$Variable[ph.table$Variable == "treatmnt"] <- "Treatment Modality"
ph.table$Variable[ph.table$Variable == "GLOBAL"] <- "Global Test"
ph.table 


#plot
plot(ph.test, main="Schoenfeld Residuals for Treatment")


#Question 9: PH test using Time Dependent Viariable
#time dependent cox model 
#we should create a numeric version of treatment again
esrd$treatmnt_num=as.numeric(esrd$treatmnt=="Haemodialysis")
tdc_model=coxph(surv_obj~treatmnt+tt(treatmnt_num), data=esrd,
                                     tt=function(x,t,...)
                                       {x*log(t+0.001)}) #added 0.001 so we dont have any potential zeroes inside the logarithm due to calculations

tdc_sum=summary(tdc_model)
tdc_sum
# Get Coefficients and P-values
tdc_sum <- summary(tdc_model)
coefs <- tdc_sum$coefficients[, "coef"]
pvals <- tdc_sum$coefficients[, "Pr(>|z|)"]
# Get Confidence Intervals for the COEFFICIENTS (not HRs)
cis <- confint(tdc_model)

#  Create the Table
tdc_df <- data.frame(
  Variable = c("Haemodialysis", "Haemodialysis: log(Time)"),
  `Coefficient` = sprintf("%.3f", coefs[c(1, 2)]),
  `95% CI` = sprintf("[%.3f, %.3f]", cis[c(1, 2), 1], cis[c(1, 2), 2]),
  `P-value` = format.pval(pvals[c(1, 2)], eps = 0.001)
)
tdc_df

####Question 10: Addressing PH Violation via Piecewise Cox Model####
#We choose as cut point the 1 year based on the LML plot and the early drop we've seen
cut_point=1

#we will include in the model treatmnt (base effect for 0-1 year)
#will include tt(treatmnt) as additional effect for 1+ year
cox_piecewise= coxph(surv_obj~treatmnt_num+tt(treatmnt_num),data = esrd,
                     tt=function(x,t,...){x*(t>1)}) #The function x*(t>1->(0,1} and (1,inf))) works like the indicator variable X(t) in Collett
model_sum_piecewise=summary(cox_piecewise)
model_sum_piecewise
#GET EARLY EFFECT 
#Interval 1 (0-1 Year) 
#This is just the main effect 'treatmnt_num'
hr_early <- model_sum_piecewise$conf.int[1, "exp(coef)"]
ci_early_lower <- model_sum_piecewise$conf.int[1, "lower .95"]
ci_early_upper <- model_sum_piecewise$conf.int[1, "upper .95"]

#CALCULATE LATE EFFECT (Manual Calculation Required)
#This is (Main Effect + Interaction Effect)
coefs <- coef(cox_piecewise)
vcov_mat <- vcov(cox_piecewise)

# We sum the coefficients
beta_late <- coefs[1] + coefs[2]

# We calculate SE for the sum: sqrt(Var1 + Var2 + 2*Cov12)
var_late <- vcov_mat[1, 1] + vcov_mat[2, 2] + 2 * vcov_mat[1, 2]
se_late <- sqrt(var_late)

hr_late <- exp(beta_late)
ci_late_lower <- exp(beta_late - 1.96 * se_late)
ci_late_upper <- exp(beta_late + 1.96 * se_late)

# Results Table
piecewise_results_df <- data.frame(
  `Time Interval` = c("0-1 Year", "1+ Years"),
  `Hazard Ratio` = sprintf("%.3f", c(hr_early, hr_late)),
  `95% CI` = c(sprintf("[%.3f, %.3f]", ci_early_lower, ci_early_upper),
               sprintf("[%.3f, %.3f]", ci_late_lower, ci_late_upper)),
  `Effectiveness (RRR)` = sprintf("%.1f%%", c((1 - hr_early)*100, (1 - hr_late)*100))
)
piecewise_results_df

#Question 11: Multivariable Cox Model
#We start with the main effects model adjusting for all baseline covariates found in Question 1

cox_multi=coxph(surv_obj~treatmnt+age+sex+copd+dm+hypert+heart+liver+neoplasia+vascular,data = esrd)
sum_cox_multi=summary(cox_multi)
sum_cox_multi


#goodness of fit with cox-snell resids
#We calculate the Nelson-Alen estimator of the cumulative hazard of the residuals

esrd$dth_num=as.numeric(esrd$dth=="Died") #we need dth as numeric again
#Calculate the martingale resids
martingale_res=residuals(cox_multi,type="martingale")
#calc cox-snell resids: cs= delta- martingale
cox_snell=esrd$dth_num-martingale_res
#Nelson Alen estimator to the resids
fit_cs=survfit(Surv(cox_snell,esrd$dth_num)~1)
#plot H(t) vs t
plot(fit_cs$time, -log(fit_cs$surv), type='l', lwd=2,
     xlab="Cox-Snell Residual", 
     ylab="Estimated Cumulative Hazard",
     main="Cox-Snell Residuals (Model Fit Check)")
abline(0, 1, col="red", lty=2) # The model fits if data follows this line


#interactions

cox_interactions <- coxph(surv_obj ~ treatmnt * (age + sex + copd + dm + hypert + heart + liver + neoplasia + vascular),data = esrd)

int_summary=summary(cox_interactions)$coefficients
int_confint <- summary(cox_interactions)$conf.int
#filter for only interaction terms (rows with `:`)
int_rows=grep(":", rownames(int_summary))
int_data=int_summary[int_rows,]
int_ci <- int_confint[int_rows, ]

#create data frame
int_table <- data.frame(
  Raw_Name = rownames(int_data),
  `Hazard_Ratio` = sprintf("%.3f", int_data[, "exp(coef)"]),
  `CI_95` = sprintf("[%.3f, %.3f]", int_ci[, "lower .95"], int_ci[, "upper .95"]),
  `P_value` = format.pval(int_data[, "Pr(>|z|)"], eps = 0.001)
)
# Clean row names for the report
int_table$Interaction <- int_table$Raw_Name
int_table$Interaction <- gsub("treatmntHaemodialysis:", "Treatment x ", int_table$Interaction)
int_table$Interaction <- gsub("Yes", "", int_table$Interaction)
int_table$Interaction <- gsub("Male", "", int_table$Interaction)

# Select final columns
final_int_table <- int_table[, c("Interaction", "Hazard_Ratio", "CI_95", "P_value")]
colnames(final_int_table) <- c("Interaction", "Hazard Ratio", "95% CI", "p-value")
final_int_table

# Re-assess Proportional Hazards (PH)
ph_test_multi <- cox.zph(cox_multi)
ph_test_multi_res=ph_test_multi$table

ph_multi_table=data.frame(
  Covariate= rownames(ph_test_multi_res),
  `Chi-squared`=ph_test_multi_res[,"chisq"],
  `P-value`=format.pval(ph_test_multi_res[,"p"],eps=0.001)
)
ph_multi_table$Covariate[ph_multi_table$Covariate == "treatmnt"] <- "Treatment"
ph_multi_table$Covariate[ph_multi_table$Covariate == "GLOBAL"] <- "Global Test"
ph_multi_table

# Graphical PH Check (Schoenfeld)
plot(ph_test_multi[1], #treatmnt
     main = "Schoenfeld Residuals for Treatment (Adjusted Model)",
     xlab = "Time", ylab = "Beta(t) for Treatment")
#  Final Model: Adjusted Piecewise Cox Model
# We use 'treatmnt_num' for the time interaction
cox_final <- coxph(surv_obj ~ treatmnt_num + tt(treatmnt_num) + 
                     age + sex + copd + dm + hypert + 
                     heart + liver + neoplasia + vascular, 
                   data = esrd,
                   tt = function(x, t, ...) { x * (t > 1) })


full_mod.sum=summary(cox_final)
#extract the final adjusted HRs for the table
coefs_final=coef(cox_final)
vcov_final=vcov(cox_final)

#Interval 0-1 Year
hr_early.final<- full_mod.sum$conf.int[1, "exp(coef)"]
ci_early_lower.final <- full_mod.sum$conf.int[1, "lower .95"]
ci_early.upper.final <- full_mod.sum$conf.int[1, "upper .95"]

#interval 1+ years
beta_late.final=coefs_final["treatmnt_num"]+coefs_final["tt(treatmnt_num)"]

#var=var(A) + var(B)+2Cov(A,B)
var_late_final=vcov_final["treatmnt_num","treatmnt_num"]+ vcov_final["tt(treatmnt_num)","tt(treatmnt_num)"]+
  2*vcov_final["treatmnt_num","tt(treatmnt_num)"]
se_late.final= sqrt(var_late_final)
hr_late.final=exp(beta_late.final)
ci_late_lo.final=exp(beta_late.final-1.96*se_late.final)
ci_late_hi.final=exp(beta_late.final+1.96*se_late.final)
#create the final table
final_res <- data.frame(
  `Time Interval` = c("0-1 Year", "1+ Years"),
  `Adjusted HR` = sprintf("%.3f", c(hr_early.final, hr_late.final)),
  `95% CI` = c(sprintf("[%.3f, %.3f]", ci_early_lower.final, ci_early.upper.final),
               sprintf("[%.3f, %.3f]", ci_late_lo.final, ci_late_hi.final)),
  `Adj. Effectiveness (RRR)` = sprintf("%.1f%%", c((1 - hr_early.final)*100, (1 - hr_late.final)*100))
)
final_res

#Question 12: Weibull parametric modeling
#we recall from Question 4 the km_fit model which is survfit(surv_obj ~ treatmnt, data = esrd)
#and now we fit a Weibull model
fitWei=survreg(surv_obj~treatmnt,data = esrd,dist = "weibull")
summary(fitWei)
#We calculate the shape parameter.
#in R scale (sigma)=1/Shape (p)
sigma_est=fitWei$scale
#SE for Log(scale)
sigma_se=summary(fitWei)$table["Log(scale)","Std. Error"]

#We now calculate the CI for sigma
sigma_ci_lo=exp(log(sigma_est)-1.96*sigma_se)
sigma_ci_hi=exp(log(sigma_est)+1.96*sigma_se)

#convert to shape (p=1-sigma)
#Caution: we swap the lower/upper bounds because 1/x flips the order
shape_mle=1/sigma_est
shape_ci_lo=1/sigma_ci_hi
shape_ci_hi=1/sigma_ci_lo

#Format nicely
shape_str <- sprintf("%.3f [%.3f, %.3f]", shape_mle, shape_ci_lo, shape_ci_hi)

# Now we Calculate SCALE Parameter (alpha) for Each Group 
#Alpha = exp(Linear Predictor)
#We use predict() to get the Linear Predictor and its Standard Error correctly

#Conservative Group
pred_cons <- predict(fitWei, newdata = data.frame(treatmnt = "Conservative"), type = "lp", se.fit = TRUE)
scale_cons_est <- exp(pred_cons$fit[1])
scale_cons_lo <- exp(pred_cons$fit[1] - 1.96 * pred_cons$se.fit[1])
scale_cons_hi <- exp(pred_cons$fit[1] + 1.96 * pred_cons$se.fit[1])
scale_cons_str <- sprintf("%.3f [%.3f, %.3f]", scale_cons_est, scale_cons_lo, scale_cons_hi)

#Haemodialysis Group
pred_haem <- predict(fitWei, newdata = data.frame(treatmnt = "Haemodialysis"), type = "lp", se.fit = TRUE)
scale_haem_est <- exp(pred_haem$fit[1])
scale_haem_lo <- exp(pred_haem$fit[1] - 1.96 * pred_haem$se.fit[1])
scale_haem_hi <- exp(pred_haem$fit[1] + 1.96 * pred_haem$se.fit[1])
scale_haem_str <- sprintf("%.3f [%.3f, %.3f]", scale_haem_est, scale_haem_lo, scale_haem_hi)


#Create and Print the Table 
mle_table <- data.frame(
  Parameter = c("Shape (p)", "Scale (alpha) - Cons.", "Scale (alpha) - Haem."),
  `MLE Estimate (95% CI)` = c(shape_str, scale_cons_str, scale_haem_str)
)
mle_table

#transform coefs
beta_ph=-coef(fitWei)/fitWei$scale
beta_ph

#Calc parameters (shape & scale)

shape_mle=1/fitWei$scale
scale_cons=exp(coef(fitWei)["(Intercept)"])
scale_haem=exp(coef(fitWei)["(Intercept)"]+coef(fitWei)["treatmntHaemodialysis"])
scale_cons
scale_haem
#plot Weibull vs KM
plot(km_fit,mark.time = F,
     main="Predicted survival for Weibull model vs KM",
     xlab = "Time (years)", ylab = "Survival Probability",lty = 1:2,col = c("blue","red"),
     lwd = 2)
#Add curves fitted by the Weibull model
pct= seq(0,1,by=0.001)
lines(predict(fitWei,newdata = data.frame(treatmnt="Conservative"),type = "quantile",
              p=pct),1-pct,lty=3,col="green",lwd=2)
lines(predict(fitWei,newdata = data.frame(treatmnt="Haemodialysis"),type = "quantile",
              p=pct),1-pct,lty=4,col="orange")
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("KM: Conservative","KM: Haemodialysis","Weibull: Conservative","Weibull: Haemodialysis"),ncol = 2)
#lml plot again: we are not going to re add it to the report as it's already done in Q8
plot(km_fit,mark.time = F,main = "Estimated log cumulative hazard vs log(time)", 
     xlab = "Time (years on log-scale)",ylab = "Log cumulative hazard",lty = 1:2,
     col = c("blue","red"), fun = "cloglog")

#12 (b) Multivariable Weibull Model

#We fit a Multivariable Weibull
#and we use the same covariates as in Question 11
wei_multi <- survreg(surv_obj ~ treatmnt + age + sex + copd + dm + 
                       hypert + heart + liver + neoplasia + vascular, 
                     data = esrd, dist = "weibull")
summary(wei_multi)
#Create "Average Patient" Data for Prediction
#We hold all covariates constant at their mean/mode to isolate the treatment effect
new_data <- data.frame(
  treatmnt = c("Conservative", "Haemodialysis"),
  age = mean(esrd$age),
  sex = "Male",       # Most common sex
  copd = "No",        # Most common
  dm = "No",
  hypert = "Yes",     # Most common
  heart = "Yes",
  liver = "No",
  neoplasia = "No",
  vascular = "No"
)
#Predict Survival Probabilities
#We predict time points (quantiles) for survival probabilities 0.99 down to 0.01
probs <- seq(0.01, 0.99, by = 0.01)
pred_times <- predict(wei_multi, newdata = new_data, type = "quantile", p = probs)
surv_probs <- 1 - probs # Convert CDF to Survival (S(t))

#Plot the Predicted Curves
plot(pred_times[2, ], surv_probs, type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "Time (Years)", ylab = "Predicted Survival Probability",
     main = "Predicted Survival (Multivariable Weibull Model)",
     xlim = c(0, 10)) # Limit to 10 years for readability
lines(pred_times[1, ], surv_probs, col = "blue", lwd = 2, lty = 1)

legend("topright", legend = c("Conservative", "Haemodialysis"),
       col = c("blue", "red"), lty = 1, lwd = 2, bty = "n")
