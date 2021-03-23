args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(RColorBrewer)


name_input  = args[1]
name_output = args[2]
print(name_input)


data <- read.csv(name_input, sep="\t")
head(data)


data$y <- factor(data$y, levels = data$y)
p <- ggplot(data, aes(x, y, fill=z, height=1), width=50) + 
  geom_tile() +
  theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face="plain", color="#222222", 
                           size=12, angle=0))+
  coord_fixed() +
  scale_fill_gradientn( colours = brewer.pal(n= 9, name="YlOrRd")) +

  labs(y='TRANSPOSABLE ELEMENT', x='', fill='INSERTION') +
  theme(text=element_text(size=14), plot.title = element_text(size=30, hjust = 0.5)) 
p
ggsave(name_output)

