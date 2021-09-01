args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(RColorBrewer)

name_input = args[1]
name_output = args[2]
print(name_input)


data <- read.csv(name_input, sep="\t")
head(data)

#reordonner 
piclust = c('pld-:-jing', 'dip1-:-178783', 'kua-:-spir', 'scro-:-nd-aggg', 'nrm-:-ago3', 'ago3-:-nvd', 'cg41099-:-aux', 'gprk1-:-ir41a', 'su(f)-:-395000', 'onecut-:-unc-13', 'fog-:-fuctc', 'myo81f-:-myo81f', 'cg17683-:-gprk1', 'cht10-:-cg12567', 'l(2)41ab-:-cg17691', 'unc-13-:-camkii', 'rpl5-:-cg40006', 'mfs17-:-cg41378', 'cadps-:-70000', 'clamp-:-marf1', 'cr45227-:-rya', 'dip-lambda-:-dip-lambda', 'sxc-:-znt41f', 'cg41378-:-scp1', 'lovit-:-fasn3', 'cg12567-:-tim23', 'cg9380-:-55000', 'his-psi:cr33867-:-eef2', 'fasn3-:-cr41320', 'eif4b-:-65000', 'zfh2-:-gat', 'kto-:-su(z)12', 'dip-lambda-:-55000', 'cr40190-:-cg40191', 'cg40006-:-cg40006', 'rangap-:-cg10194', 'rdga-:-rdga', 'cg33552-:-vhasfd', 'cg8407-:-rps11', 'smr-:-cg32647', 'dbp80-:-dbp80', 'cyt-c1l-:-cg11951')
data$y <- factor(data$y, levels = rev(piclust))

#couleur gradiant blue
ggplot(data, aes(x, y, fill=z, width=0.95, height=0.95)) + 
  geom_tile()+

  scale_fill_gradientn( colours = brewer.pal(n= 9, name="GnBu"))+

  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 55, vjust = 1, size = 12, hjust = 1))+
  coord_fixed()+

  labs(y='PiCluster', x='ET', fill='Insertion')+
  theme(text=element_text(size=12)) 

#sp+scale_color_gradientn(colours = rainbow(19))
ggsave(name_output)
