library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(phylosignal)
library(phangorn)


# ÁRVORE CONSENSO

MCC<-maxCladeCred(sapajus.tree,rooted=TRUE)

## JUNTANDO OS DADOS COM A FILOGENIA
# LIGHT
# criando uma matriz com os dados que tu vai usar pra representar na filogenia (no caso é a coluna 1 [,1])
Size <- as.matrix(pgls_data)[,4]

# estimativa do estado ancestral
fit <- phytools::fastAnc(MCC, Size, vars=TRUE, CI=TRUE)

# juntando os dados
td <- data.frame(node = nodeid(MCC, names(Size)), trait = Size)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)

tree <- full_join(MCC, d, by = 'node')

??fastAnc

###RED
Redness <- as.matrix(pgls_data)[,2]

fit <- phytools::fastAnc(MCC, Redness, vars=TRUE, CI=TRUE)

td <- data.frame(node = nodeid(MCC, names(Redness)), trait = Redness)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)

tree <- full_join(MCC, d, by = 'node')

# FILOGENIA PARA O LIGHTNESS
windows()

ggtree(tree, aes(color = trait)) +
  geom_tree(size = 1.5) +
  scale_colour_viridis_c(option = "viridis") +  # Set discrete viridis palette
  theme_minimal() +
  xlab("Branch Length") +
  geom_tiplab(aes(color = trait))  # Optional tip labels (redundant with tip_labels)

# FILOGENIA PARA O REDNESS# FILOGENIA PARApgls_data O REDNESS
ggtree(tree, layout = 'circular', size = 2.5) + geom_tree(aes(color = trait), continuous = 'colour', size = 2) + scale_colour_continuous(low = "gainsboro", high = "darkred", name = "Redness") + geom_tiplab(aes(color = trait), hjust = -.1, size = 1.8) + theme(legend.position = c(.01, .90))

# variação ao longo do tempo
#-- light
ggtree(tree, aes(color=trait), continuous = 'colour', yscale = "trait", size = 1.7) + 
  scale_color_continuous(low = "black", high = "gainsboro", name = "Lightness")+ theme_minimal()

#-- red
ggtree(tree, aes(color=trait), continuous = 'colour', yscale = "trait", size = 1.7) + 
  scale_color_continuous(low = "gainsboro", high = "darkred", name = "Redness")+ theme_minimal()


####---- trait map usando viridis 
# viridis pro red
ggtree(tree, aes(color=trait), continuous = 'colour', yscale = "trait", size = 1.7) + 
  scale_color_viridis_c(option = "inferno", name = "Redness")+ theme_minimal()


#viridis pro light
ggtree(tree, aes(color=trait), continuous = 'colour', yscale = "trait", size = 1.7) + 
  scale_color_viridis_c(option = "magma", name = "Lightness")+ theme_minimal()




