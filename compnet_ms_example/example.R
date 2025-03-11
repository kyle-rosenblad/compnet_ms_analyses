# Load and attach package
library(compnet) 

# Prep data, use more aesthetically pleasing version of trait name
mytraits <- ex_traits[c("ndtrait")]
names(mytraits) <- c("ND Trait")

# Build model:
mod1 <- buildcompnet(presabs = ex_presabs, spvars_multi_int = mytraits, rank=1)

# Plot results
library(patchwork)
example_violin <- fixedeff_violins(mod=mod1)
example_scatter <- scatter_interaction(mod = mod1, xvar = "ND Trait", ymax=0.4)
example_violin + example_scatter
library(ggplot2)
ggsave("example_plots.png", height=4, width=8)
