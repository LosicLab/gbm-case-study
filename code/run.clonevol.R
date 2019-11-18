# ---------------------------------------------------
# Title: Clonal Evolution Analysis
#        for GBM Case Study (Part III)
# Author: Paula Restrepo, paularstrpo@gmail.com
# Last updated: 11/04/2019
# ---------------------------------------------------

# ---------------------------------------------------
# 1) Run sciclone
# ---------------------------------------------------

libs <- c("tidyverse", "clonevol")
lapply(libs, library, character.only=TRUE)
rm(libs)
load('data/sciclone.results.rdata')
# ---------------------------------------------------
# 1) Run Clonevol
# ---------------------------------------------------


## prepare clonevol input
vafs <- sc.variants %>%
        dplyr::select(`cluster`, `X10647A.vaf`, `X10647B.vaf`, `X10647C.vaf`, `X9260B.vaf`) %>%
        mutate(`Recurrent_A`=`X10647A.vaf`,
               `Recurrent_B`=`X10647B.vaf`,
               `Recurrent_C`=`X10647C.vaf`,
               `Primary`=`X9260B.vaf`) %>%
    dplyr::select(- `X10647A.vaf`, - `X10647B.vaf`,
               -`X10647C.vaf`, -`X9260B.vaf`)

vafs <- vafs[!is.na(vafs$cluster) & vafs$cluster > 0,]

# prepare sample grouping
sample.groups <- as.factor(c('P','R', 'R', 'R'))
names <- colnames(vafs[,2:5])
names <- names[c(4,1:3)]
names(sample.groups) <- names
clone.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e')

x <- vafs
pdf('box.pdf', width = 6, height = 16)
pp <- plot.variant.clusters(x,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'black',
                            vaf.col.names = names,
                            vaf.limits = 70,
                            sample.title.size = 20,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.shape = 1,
                            jitter.color = c('#f8e356', '#fe9536','#8d4891','#999793'),
                            jitter.size = 3,
                            jitter.alpha = 1,
                            jitter.center.method = 'median',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'none',
                            order.by.total.vaf = TRUE)
dev.off()

# plot clusters pairwise-ly
plot.pairwise(x, col.names = names,
              colors = clone.colors)

# plot mean/median of clusters across samples (cluster flow)
pdf('cluster_flow.pdf', width=12, height=10)
plot.cluster.flow(x, vaf.col.names = names,
                  colors = clone.colors)
dev.off()

# infer consensus clonal evolution trees
y = infer.clonal.models(variants = x,
                        cluster.col.name = 'cluster',
                        vaf.col.names = names,
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 10000,
                        clone.colors = clone.colors,
                        founding.cluster = 1,
                        cluster.center = 'median',
                        min.cluster.vaf = 0,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05)


# prepare branch-based trees
clonal_models <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

# plot variant clusters, bell plots, cell populations, and trees
plot.clonal.models(clonal_models,
                   models.to.plot = 1,
                   # box plot parameters
                   # box.plot = TRUE,
                   # fancy.boxplot = TRUE,
                   # fancy.variant.boxplot.jitter.alpha = 1,
                   # fancy.variant.boxplot.jitter.center.color = 'grey50',
                   # fancy.variant.boxplot.base_size = 12,
                   # fancy.variant.boxplot.plot.margin = 1,
                   # fancy.variant.boxplot.order.by.total.vaf = TRUE,
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'black',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   nonzero.cell.frac.clone.border.color = 'black',
                   zero.cell.frac.clone.border.color = 'black',
                   cell.frac.side.arrow.width=1.25,
                   #  node-based consensus tree parameters
                   # merged.tree.plot = TRUE,
                   # tree.node.label.split.character = '\n',
                   # tree.node.num.samples.per.line=1,
                   # tree.node.shape = 'circle',
                   # tree.node.size = 30,
                   # merged.tree.node.size.scale = 1.25,
                   # merged.tree.node.text.size.scale = 2.5,
                   # merged.tree.cell.frac.ci = FALSE,
                   # # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = FALSE,
                   # # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.1,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   #output figure parameters
                   out.dir = 'figures/',
                   out.prefix = 'clonal_architecture_final',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 8,
                   panel.widths = c(1,1,1.5))
