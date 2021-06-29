##graph tune results

#load library
library(ggplot2)

#qualified results
#load data
tune.results <- read.csv("./res/tune.results.csv", stringsAsFactors = FALSE, header = TRUE)

#create ggplot
graph <- ggplot() +
  geom_line(data = tune.results, aes(x = ntree.per.trial, y = average.n.unique.elem.per.param, color = "black"), linetype = "solid", size = 1) +
  geom_line(data = tune.results, aes(x = ntree.per.trial, y = average.hamming.dist.per.param, color = "blue"), linetype = "solid", size = 1) +
  geom_ribbon(data = tune.results, aes(x = ntree.per.trial, ymin = average.n.unique.elem.per.param - sd.n.unique.elem.per.param, ymax = average.n.unique.elem.per.param + sd.n.unique.elem.per.param), fill = "grey") + 
  geom_ribbon(data = tune.results, aes(x = ntree.per.trial, ymin = average.hamming.dist.per.param - sd.hamming.dist.per.param, ymax = average.hamming.dist.per.param + sd.hamming.dist.per.param), fill = "blue") +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 345000) +
  labs(x = "ntree", y = "", color = "") +
  scale_color_manual(labels = c("Average Number of Pairwise Unique Elements", "Average Pairwise Hamming Distance"), values = c("black", "blue")) +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top", axis.text = element_text(color = "black"))

graph
ggsave("./fig/tune.results.eps", plot = graph)

#refined results
#load data
tune.results.refined <- read.csv("./res/tune.results.refined.csv", header = TRUE)

graph.refined <- ggplot() + 
  geom_line(data = tune.results.refined, aes(x = ntree, y = mean), linetype = "solid", size = 1) +
  labs(x = "ntree", y = "Estimated Out-of-bag Error Rate", color = "black") + 
  scale_color_manual(labels = c("Estimated Out-of-bag Error Rate"), values = c("black")) +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top", axis.text = element_text(color = "black"))

graph.refined
ggsave("./fig/refined.tune.results.eps", plot = graph.refined)
