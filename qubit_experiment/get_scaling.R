library(data.table)
library(ggplot2)

get_LMC_control <- function(filename, ctl_type='LMC'){
  df <- data.frame(fread(filename))
  df$dilution_factor <- (df$water_input_ul+df$Stock_input_ul)/df$Stock_input_ul
  df$mean_qubit <- rowMeans(df[, grepl("Qubit", colnames(df))])
  df$estimated_stock_concentration <- df$dilution_factor * df$mean_qubit
  df$control <- 'HMC'
  df$control[grepl("L", df$Sample)] <- 'LMC'
  #ggplot(data=df) + geom_point(aes(x=dilution_factor, y= estimated_stock_concentration, color = control), size=3)+theme_bw()
  LMC <- mean(df[(grepl(ctl_type, df$control) & df$mean_qubit<9 & df$mean_qubit >1 &!is.na(df$mean_qubit)), ]$estimated_stock_concentration)
  return(LMC)
}

LMC1 <- get_LMC_control('qubit_values_1.txt')
LMC2 <- get_LMC_control('qubit_values_2.txt')
HMC2 <- get_LMC_control('qubit_values_2.txt', ctl_type = 'HMC')

saveRDS(LMC1, file = 'LMC_true_concentration_1.rds')
saveRDS(LMC2, file = 'LMC_true_concentration_2.rds')
saveRDS(HMC2, file = 'HMC_true_concentration_2.rds')
# LMC is always 8ul of stock, diluted to 2000ul total
#LMC_used_value <- LMC_used_value*8/2000
#saveRDS(LMC_used_value, file = 'LMC_true_concentration.rds')
