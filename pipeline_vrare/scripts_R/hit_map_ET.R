args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(RColorBrewer)
#library(extrafont) 
#font_import() 
#loadfonts(device = "postscript")

name_file = args[1]
print(name_file)

data <- read.csv(name_file, sep="\t")
head(data)

data$y <- factor(data$y, levels = data$y)

spl1  = unlist(strsplit(name_file, "/"))
spl2  = unlist(strsplit(spl1[1], "\\."))
spl   = unlist(strsplit(spl2[1], "_"))

print(spl[1])

ggplot(data, aes(x, y, fill=z, height=1), width=50) + 
  geom_tile()+

  theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face="plain", color="#222222", 
                           size=12, angle=0))+
  coord_fixed()+

  scale_fill_gradientn(breaks = c(0, 25, 50, 75, 93), colours = brewer.pal(n= 9, name="YlOrRd"))+

  labs(y='Transposables Element', x='', fill='NEW INSERTION')+
  theme(text=element_text(size=14), plot.title = element_text(size=30, hjust = 0.5)) 

ggsave(paste(spl[1], "_hit_map.png", sep=""))

