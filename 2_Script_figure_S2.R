##### FIGURE S2: read length distribution comparison between WT, UPOP and ePRINT

#to run in R
data <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig2/peak_lengths_summary.txt", quote="\"", comment.char="")


#normalization for group
data$V1[data$V3 == "Cell"] <- (data$V1[data$V3 == "Cell"]/sum(data$V1[data$V3 == "Cell"]))  *1000000
data$V1[data$V3 == "CM"] <- (data$V1[data$V3 == "CM"]/sum(data$V1[data$V3 == "CM"]))  *1000000
data$V1[data$V3 == "ePRINT"] <- (data$V1[data$V3 == "ePRINT"]/sum(data$V1[data$V3 == "ePRINT"]))  *1000000
data$V1[data$V3 == "UPOP-seq"] <- (data$V1[data$V3 == "UPOP-seq"]/sum(data$V1[data$V3 == "UPOP-seq"]))  *1000000


ggplot(data, aes(x = V2, y = V1 , color = V3, fill = V3)) +
  geom_area(alpha = c(0.4), position = 'identity', size = 1.2)+
  scale_fill_manual(values = c("whitesmoke", "whitesmoke","#8c959d", "#babfc4")) +
  scale_color_manual(values = c("#e6a3c3", "#62c3ef","#8c959d", "#babfc4")) +  # Assigning fill colors
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
#read_length_comparison
