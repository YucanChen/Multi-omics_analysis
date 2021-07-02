library(openxlsx)
library(stringr)
library(ggplot2)
library(dplyr)

# Using Only BP&REACTOME term for g:Profiler enrichment
DEmiRNA_regulated_DOWN <- read.csv("gProfiler_DEmiRNA_regulated_DOWN_GOBPReactome.csv")
DEmiRNA_regulated <- read.csv("gProfiler_DEmiRNA_regulated_GOBPReactome.csv")
DEmiRNA_regulated_UP <- read.csv("gProfiler_DEmiRNA_regulated_UP_GOBPReactome.csv")
H3K36me3_regulated <- read.csv("gProfiler_H3K36me3_regulated_GOBPReactome.csv")
peakcoverDEmiRNA_regulated_DOWN <- read.csv("gProfiler_peakcoverDEmiRNA_regulated_DOWN_GOBPReactome.csv")
peakcoverDEmiRNA_regulated <- read.csv("gProfiler_peakcoverDEmiRNA_regulated_GOBPReactome.csv")
### DEmiRNA/H3K36me3/peakcoverDEmiRNA
# add group marks
DEmiRNA_regulated_DOWN$group <- rep("DEmiRNA-regulated down",5)
DEmiRNA_regulated_UP$group <- rep("DEmiRNA-regulated up",4)
DEmiRNA_regulated$group <- rep("DEmiRNA-regulated all",2)
H3K36me3_regulated$group <- rep("H3K36me3-regulated",5)
peakcoverDEmiRNA_regulated$group <- rep("peakcovered DEmiRNA-regulated all",1)
peakcoverDEmiRNA_regulated_DOWN$group <- rep("peakcovered DEmiRNA-regulated down",117)

# combine BP term (intersection>=2, catagory:5-350)
combined_BP_plot <- H3K36me3_regulated[which((H3K36me3_regulated$source=="GO:BP") & (H3K36me3_regulated$intersection_size >= 2)),]
tmp <- DEmiRNA_regulated[which((DEmiRNA_regulated$source=="GO:BP") & (DEmiRNA_regulated$intersection_size >= 2)),]
combined_BP_plot <- bind_rows(combined_BP_plot,tmp)
tmp <- DEmiRNA_regulated_DOWN[which((DEmiRNA_regulated_DOWN$source=="GO:BP") & (DEmiRNA_regulated_DOWN$intersection_size >= 2)),]
combined_BP_plot <- bind_rows(combined_BP_plot,tmp)
tmp <- DEmiRNA_regulated_UP[which((DEmiRNA_regulated_UP$source=="GO:BP") & (DEmiRNA_regulated_UP$intersection_size >= 2)),]
combined_BP_plot <- bind_rows(combined_BP_plot,tmp)
tmp <- peakcoverDEmiRNA_regulated[which((peakcoverDEmiRNA_regulated$source=="GO:BP") & (peakcoverDEmiRNA_regulated$intersection_size >= 2)),]
combined_BP_plot <- bind_rows(combined_BP_plot,tmp)
tmp <- peakcoverDEmiRNA_regulated_DOWN[which((peakcoverDEmiRNA_regulated_DOWN$source=="GO:BP") & (peakcoverDEmiRNA_regulated_DOWN$intersection_size >= 2)),]
combined_BP_plot <- bind_rows(combined_BP_plot,tmp)

## Draw the plot
ggplot(data=combined_BP_plot)+
  geom_bar(aes(x=term_name,y=negative_log10_of_adjusted_p_value,fill=group),stat='identity',position = "dodge",width=0.7) + 
  scale_fill_manual(values = c("#CC33CC","#6666CC","#660099","#996600","#669933"))+ 
  coord_flip() +
  scale_x_discrete(labels=function(x) str_wrap(x, width=30)) +
  xlab("") +
  ylab("-log10(adjusted P-value)") +
  theme(
    axis.text.x=element_text(color="black",size=rel(1.1)),
    axis.text.y=element_text(color="black", size=rel(1.1)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    legend.text=element_text(color="black",size=rel(0.9)),
    legend.title = element_text(color="black",size=rel(1.0))
  ) +
  theme(panel.background = element_rect(color = 'grey60', fill = 'transparent'))+
  geom_hline(yintercept = -1*log10(0.05),colour = "#333399",lty = 2)

ggplot(data=combined_REAC_plot)+
  geom_bar(aes(x=term_name,y=negative_log10_of_adjusted_p_value,fill=group),stat='identity',position = "dodge",width=0.6) + #position_dodge(0.9)
  scale_fill_manual(values = c("#CC33CC","#CC0033","#660099"))+
  coord_flip() +
  scale_x_discrete(labels=function(x) str_wrap(x, width=30)) +
  xlab("") +
  ylab("-log10(adjusted P-value)") +
  theme(
    axis.text.x=element_text(color="black",size=rel(1.1)),
    axis.text.y=element_text(color="black", size=rel(1.1)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    legend.text=element_text(color="black",size=rel(0.9)),
    legend.title = element_text(color="black",size=rel(1.0))
  ) +
  geom_hline(yintercept = -1*log10(0.05),colour = "#333399",lty = 2)+
  theme(panel.background = element_rect(color = 'grey60', fill = 'transparent'))
