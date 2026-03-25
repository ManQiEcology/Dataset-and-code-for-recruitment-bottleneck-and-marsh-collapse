library(here)
library(readxl)
library(caret)

setwd(here())



###################################################
#DI_2018
Confs_DI_2018<-read_excel("Data/Vegetation classification/Deal Island/2018/DI_2018_accuracty table.xlsx")

combined_levels <- union(levels(factor(Confs_DI_2018$Ground_truth_IDENT)), levels(factor(Confs_DI_2018$Classified_IDENT)))
Confs_DI_2018$Ground_truth_IDENT <-factor(Confs_DI_2018$Ground_truth_IDENT, levels = combined_levels)
Confs_DI_2018$Classified_IDENT<-factor(Confs_DI_2018$Classified_IDENT, levels = combined_levels)
Confs_DI_2018$Post_processed_IDENT<-factor(Confs_DI_2018$Post_processed_IDENT, levels = combined_levels)
conf_matrix <- confusionMatrix(Confs_DI_2018$Classified_IDENT, Confs_DI_2018$Ground_truth_IDENT)
conf_matrix_post_processing <- confusionMatrix(Confs_DI_2018$Post_processed_IDENT, Confs_DI_2018$Ground_truth_IDENT)

print(conf_matrix)
print(conf_matrix_post_processing)

# Extract the confusion matrix table
conf_matrix_table <- as.data.frame(conf_matrix$table)
conf_matrix_table <- as.data.frame(conf_matrix_post_processing$table)

# Plot the confusion matrix using ggplot2

tiff("Confusion table_DI_2018.tiff", unit="in",width = 10, height =8, res= 600,pointsize = 10)
ggplot(data = conf_matrix_table, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_text(aes(label = Freq), vjust = 1) +
  labs(x = "Actual classification", y = "Estimated classification") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()