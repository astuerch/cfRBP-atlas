
### FIGURE S5A: correlation between RNA fragment size (footprint) and RBP size (InterPro)
### FIGURE S5B: correlation between RNA fragment size (footprint) and RNA-binding domain size (InterPro)

#correlation with RBP and RBD sizes
reads <- read.table("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/3_Footprint_characterization/fragment_size/all_RBPs_master_length_summary.txt", quote="\"", comment.char="", header = T)

cor_coef1 <- cor(reads$avg_V1, reads$RBD_size) #pearson test

lm_model <- lm(avg_V1 ~ RBD_size, data = reads)
summary(lm_model)

r_squared1 <- cor_coef1^2
max_value <- max(c(max(reads$avg_V1), max(reads$RBD_size)))


ggplot(data = reads)+
  geom_point(aes(x = avg_V1, y = RBD_size), color = "black", alpha = 0.9) +
  geom_smooth(aes(x = avg_V1, y = RBD_size), method = "lm", color = "red", se = FALSE, size = 1) +
  labs(x = "Mean RNA fragment size [nt]", y = "RBD size [nt]") +
  #coord_fixed(ratio = 1) +
  #xlim(0, max_value) +  # Set x-axis limits
  #ylim(0, max_value) +  # Set y-axis limits
  theme_classic(base_size = 15, base_family = "Helvetica") +
  annotate("text", x = 16, y = 300, label = paste0("R = ", round(cor_coef1, 2)),col = "red", size = 6) +
  annotate("text", x = 16, y = 250, label = paste0("p-value = 0.97"),col = "red", size = 6)

ggplot(data = reads)+
  geom_point(aes(x = avg_V1, y = RBP_size), color = "black", alpha = 0.9) +
  geom_smooth(aes(x = avg_V1, y = RBP_size), method = "lm", color = "#489FA7", se = FALSE, size = 1) +
  labs(x = "Mean RNA fragment size [nt]", y = "RBP size [nt]") +
  #coord_fixed(ratio = 1) +
  #xlim(0, max_value) +  # Set x-axis limits
  #ylim(0, max_value) +  # Set y-axis limits
  theme_classic(base_size = 15, base_family = "Helvetica") +
  annotate("text", x = 16, y = 1990, label = paste0("R = ", round(cor_coef1, 2)),col = "#489FA7", size = 6) +
  annotate("text", x = 16, y = 1800, label = paste0("p-value = 0.2443"),col = "#489FA7", size = 6)





### FIGURE S5C: RNA fragments distribution in HEK CM and plasma

#read length from bed reliable filtered files
samtools view Plasma_detector_merged.bam | grep -E '(NH:i:1|^@)' | awk '{print \$3, length(\$10)}' | awk '{print $2, "Plasma"}' | sort | uniq -c > Plasma_detector_read_length.txt

samtools view WT_CM_smRNA_allbatches_merged.bam | grep -E '(NH:i:1|^@)' | awk '{print \$3, length(\$10)}' | awk '{print $2, "HEK_CM"}' | sort | uniq -c > WT_CM_smRNA_master_fragments_read_length.txt


cat Plasma_detector_read_length.txt WT_CM_smRNA_master_fragments_read_length.txt  > all_read_length.txt


#r
data <- read.table("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/2_Footprint_biofluids/all_read_length.txt", quote="\"", comment.char="")

#normalization for biotypes
data$V1[data$V3 == "HEK_CM"] <- (data$V1[data$V3 == "HEK_CM"]/sum(data$V1[data$V3 == "HEK_CM"]))  *1000000
data$V1[data$V3 == "Plasma"] <- (data$V1[data$V3 == "Plasma"]/sum(data$V1[data$V3 == "Plasma"]))  *1000000


ggplot(data, aes(x = V2, y = V1 , color = V3, fill = V3)) +
  geom_area(alpha = c(0.2), position = 'identity', size = 1.2)+
  #scale_fill_manual(values = c("whitesmoke", "whitesmoke", "whitesmoke")) +
  scale_fill_manual(values = c( "#64C2EE", "#E79835" )) +  # Assigning fill colors
  scale_color_manual(values = c( "#64C2EE", "#E79835" )) +  # Assigning fill colors "#1C3F5A", "#6BBFB4"
  # Assigning fill colors
  #scale_x_continuous(breaks = unique(as.numeric(factor(data$V2))), labels = levels(factor(data$V2))) +
  scale_x_continuous(n.breaks = 14) +
  labs(title = "All biotypes",
       x = "Length / nt",
       y = "RPM") +
  theme_classic(base_family = "Arial")+
  theme(axis.text = element_text(size = 8),  # Adjust font size for axis text
        axis.title = element_text(size = 15),  # Adjust font size for axis titles
        plot.title = element_text(size = 20, face = "bold"))

#6x5
#read_length_distribution_bam


