################################################################################
########## Project: GMDS Summer School                                ##########
########## Purpose: Generate Forest Plots for resutls                 ##########
########## Author: Halimu Haliduola                                   ##########
########## Date: 03Sep2024                                            ##########
########## Updates:                                                   ##########
################################################################################



#install.packages("forestplot")
library(forestplot)

base_data <- tibble::tibble(mean  = c(2.03, 2.15, 2.66, 2.36, 1.52, 3.36),
                            lower = c(1.49, 1.62, 2.18, 1.84, 0.93, 2.76),
                            upper = c(2.58, 2.68, 3.14, 2.89, 2.10, 3.96),
                            Method = c("True", "UBR", "Random Forest", "MI under MAR",
                                      "MI under MNAR", "No Imputation"),
                            Test_Grp_Mean = c("12.01", "12.14", "12.66", "12.41", "11.56", "13.28"),
                            Ctrl_Grp_Mean = c("9.98",  "9.98",  "10.00", "10.04", "10.05", "9.92"),
                            N_Total = c("600", "600", "600", "600", "600", "445"))


base_data |>
  forestplot(labeltext = c(Method, Test_Grp_Mean, Ctrl_Grp_Mean, N_Total),
             clip = c(0.5, 5.0),
             vertices= TRUE,
             xlog = FALSE) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(Method = c("", "Method"),
                Test_Grp_Mean = c("Test", "(LS Mean)"),
                Ctrl_Grp_Mean = c("Control", "(LS Mean)"),
                N_Total = c("N", "Total")) |>
  fp_set_zebra_style("#EFEFEF")



