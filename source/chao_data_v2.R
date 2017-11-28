# |----------------------------------------------|
# | Project: Curcumin Treatment of AOM-DSS Model |
# | Study ID:                                    |
# | Scientist: Yue, Renyi Wu                     |
# | Data Analysis: Renyi Wu, Davit Sargsyan      |
# | Created: 11/06/2017                          |
# |----------------------------------------------|
# Header----
# require(devtools)
# install_github("ggbiplot", "vqv")

require(data.table)
require(ggplot2)
require(ggbiplot)
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)

# Part I: Load data and annotate----
peakAnno1 <- annotatePeak(peak = "data/chao/Methyl_pten.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                          annoDb = "org.Mm.eg.db")

dt1 <- data.table(as.data.frame(peakAnno1@anno@ranges),
                  as.data.frame(peakAnno1@anno@elementMetadata@listData))
dt1

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]

# OPTIONAL: set NAs to zero, add 1 to all data points
# dt1[, 5:16] <- lapply(dt1[, 5:16],
#                       function(a) {
#                         a[is.na(a)] <- 0
#                         a <- a + 1
#                       })
dt1

# Part III: % Methylation----
# Rename all samples using treatment names
sName <- colnames(dt1)[5:28]
sName <- unique(substr(sName,
                       1,
                       nchar(sName) - 2))
sName

t1 <- data.table(dt1[, c(1:4, 29, 30, 35, 37, 39)], 
                 dt1[, seq(6, 
                           28,
                           2), 
                     with = FALSE]/ 
                   dt1[, seq(5, 
                             27,
                             2), 
                       with = FALSE])
names(t1)[10:ncol(t1)] <- sName
t1
gc()

# Part III: PCA----
# Remove all incomplete 
tmp <- t1[, 10:21,
          with = FALSE]
ndx <- is.na(rowSums(tmp))
sum(ndx)
# Excluded 104,621 regions for now
tmp <- t(tmp[!ndx, ])
dim(tmp)
colnames(tmp) <- t1$start[!ndx]

m.pca <- prcomp(tmp)
summary(m.pca)
plot(m.pca)

# NEW (11/11/2017): Keep only the most important variables (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2

nobs.factor <- sqrt(nrow(m.pca$x) - 1)
d <- m.pca$sdev
u <- m.pca$x
v <- m.pca$rotation

# Scores
df.u <- data.frame(u[, choices])
# Add grouping variable
rownames(df.u)
df.u$grp <- factor(rep(1:6, each = 2))
df.u

# Directions
df.v <- as.data.frame(v[, choices])

# Annotate
df.v$abrvName <- rownames(df.v)

# Separate top variables (largest arrows)
df.v$lgth <- sqrt(df.v$PC1^2 + df.v$PC2^2)
df.v <- df.v[order(df.v$lgth,
                   decreasing = TRUE), ]
df.v

# Top and bottom 10 genes
p0 <- ggplot(df.v[c(1:10,
                    500:510), ]) +
  geom_bar(aes(x = abrvName,
               y = lgth),
           stat = "identity") +
  scale_x_discrete("CpG Start") +
  scale_y_continuous("Axis Length") +
  theme(axis.text.x = element_text(angle = 75,
                                   hjust = 1))
p0
tiff(filename = "tmp/pca_varAxisLength.tiff",
     height = 5,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p0)
graphics.off()

# Axis labels
u.axis.labs <- paste(colnames(df.v)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m.pca$sdev[choices]^2/sum(m.pca$sdev^2)))

var.keep.ndx <- 1:10

p1 <- ggplot(data = df.v[var.keep.ndx,], 
             aes(x = PC1,
                 y = PC2)) +
  coord_equal() +
  geom_point(data = df.u,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0) +
  # geom_segment(aes(x = 0,
  #                  y = 0, 
  #                  xend = 1000*PC1,
  #                  yend = 1000*PC2),
  #              arrow = arrow(length = unit(1/2, 'picas')), 
  #              color = muted('red'),
  #              size = 1.2) +
  # geom_text(aes(label = df.v$abrvName[var.keep.ndx]),
  #           size = 5,
  #           angle = 10,
  #           hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  # scale_fill_manual(name = "Treatment",
  #                   labels = c("KO 12w",
  #                              "KO 20w",
  #                              "UA 12w",
  #                              "UA 20w",
  #                              "WT 12w",
  #                              "WT 20w"),
  #                   values = c("red",
  #                              "blue",
  #                              "green",
  #                              "black",
  #                              "white",
  #                              "yellow")) +
  ggtitle("Biplot of UA Treatmed Samples") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))
p1

tiff(filename = "tmp/pca_top_vars.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()












p1 <- ggbiplot(m1,
         obs.scale = 1,
         var.scale = 1,
         varname.size = 5,
         varname.adjust = 7,
         groups = t3$start[!ndx])

tiff(filename = "tmp/pca1.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(p1)
graphics.off()




# Part II: Compare Treatments and Batches----
# Combine first and third batches
dtc <- merge(t1, t3, by = names(t1)[1:4])

# Group averages: Batch 1----
dtc$ctrl_x <- (dtc$Control_1.x + dtc$Control_2.x)/2
dtc$aomdss_x <- (dtc$AOM_DSS_1.x + dtc$AOM_DSS_2.x)/2
dtc$cur_x <- (dtc$Cur_1.x + dtc$Cur_2.x)/2

# Fold-changes----
dtc$log2.aomdss.ctrl.x <- log2(dtc$aomdss_x/dtc$ctrl_x)
dtc$log2.aomdss.cur.x <- log2(dtc$aomdss_x/dtc$cur_x)

# Absolute differences----
dtc$diff.aomdss.ctrl.x <- dtc$aomdss_x - dtc$ctrl_x
dtc$diff.aomdss.cur.x <- dtc$aomdss_x - dtc$cur_x

# Plot proportions methylated----
hist(dtc$ctrl_x, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 1 Control")
hist(dtc$aomdss_x, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 1 AOM-DSS")
hist(dtc$cur_x, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 1 Curcumin")

# Differences in treatments within Batch 1-----
x <- dtc$ctrl_x
y <- dtc$aomdss_x
z <- dtc$cur_x

plot(y ~ x,
     xlab = "Control",
     ylab = "AOM-DSS",
     main = "Batch 1. RED:Diff <= -10% or >= 10%")
points(y[abs(dtc$diff.aomdss.ctrl.x) >= 0.1] ~ 
         x[abs(dtc$diff.aomdss.ctrl.x) >= 0.1],
       col = "red")

plot(z ~ y,
     xlab = "AOM-DSS",
     ylab = "Curcumin",
     main = "Batch 1. RED:Diff <= -10% or >= 10%")
points(z[abs(dtc$diff.aomdss.cur.x) >= 0.1] ~ 
         y[abs(dtc$diff.aomdss.cur.x) >= 0.1],
       col = "red")

# Fold-changes in treatments within Batch 1-----
x <- -log2(dtc$ctrl_x)
y <- -log2(dtc$aomdss_x)
z <- -log2(dtc$cur_x)

plot(y ~ x,
     xlab = "-log2(Control)",
     ylab = "-log2(AOM-DSS)",
     main = "Batch 1. RED:log2 <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.x) >= 1],
       col = "red")

plot(z ~ y,
     xlab = "-log2(AOM-DSS)",
     ylab = "-log2(Curcumin)",
     main = "Batch 1. RED:log2 <= -1 or >= 1")
points(z[abs(dtc$log2.aomdss.cur.x) >= 1] ~ 
         y[abs(dtc$log2.aomdss.cur.x) >= 1],
       col = "red")

# Group averages: Batch 2 Rerun----
dtc$ctrl_y <- (dtc$Control_1.y + dtc$Control_2.y)/2
dtc$aomdss_y <- (dtc$AOM_DSS_1.y + dtc$AOM_DSS_2.y)/2
dtc$cur_y <- (dtc$Cur_1.y + dtc$Cur_2.y)/2

# Differences in treatments within Batch 2 Rerun-----
# Fold-cahnges----
dtc$log2.aomdss.ctrl.y <- log2(dtc$aomdss_y/dtc$ctrl_y)
dtc$log2.aomdss.cur.y <- log2(dtc$aomdss_y/dtc$cur_y)

# Absolute differences----
dtc$diff.aomdss.ctrl.y <- dtc$aomdss_y - dtc$ctrl_y
dtc$diff.aomdss.cur.y <- dtc$aomdss_y - dtc$cur_y

# Plot proportions methylated
hist(dtc$ctrl_y, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 2 Rerun Control")
hist(dtc$aomdss_y, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 2 Rerun AOM-DSS")
hist(dtc$cur_y, 
     100, 
     xlab = "Proportion Methylated CpG",
     main = "Batch 2 Rerun Curcumin")

# Differences in treatments within Batch 2 Rerun-----
x <- dtc$ctrl_y
y <- dtc$aomdss_y
z <- dtc$cur_y

plot(y ~ x,
     xlab = "Control",
     ylab = "AOM-DSS",
     main = "Batch 2 Rerun. RED:Diff <= -10% or >= 10%")
points(y[abs(dtc$diff.aomdss.ctrl.y) >= 0.1] ~ 
         x[abs(dtc$diff.aomdss.ctrl.y) >= 0.1],
       col = "red")

plot(z ~ y,
     xlab = "AOM-DSS",
     ylab = "Curcumin",
     main = "Batch 2 Rerun. RED:Diff <= -10% or >= 10%")
points(z[abs(dtc$diff.aomdss.cur.y) >= 0.1] ~ 
         y[abs(dtc$diff.aomdss.cur.y) >= 0.1],
       col = "red")

# Fold-cahnges in treatments within Batch 2-----
x <- -log2(dtc$ctrl_y)
y <- -log2(dtc$aomdss_y)
z <- -log2(dtc$cur_y)

plot(y ~ x,
     xlab = "-log2(Control)",
     ylab = "-log2(AOM-DSS)",
     main = "Batch 2 Rerun. RED:log2 <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.y) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.y) >= 1],
       col = "red")

plot(z ~ y,
     xlab = "-log2(AOM-DSS)",
     ylab = "-log2(Curcumin)",
     main = "Batch 2 Rerun. RED:log2 <= -1 or >= 1")
points(z[abs(dtc$log2.aomdss.cur.y) >= 1] ~ 
         y[abs(dtc$log2.aomdss.cur.y) >= 1],
       col = "red")

# Fold-changes between batches----
plot(dtc$log2.aomdss.ctrl.x  ~ 
       dtc$log2.aomdss.ctrl.y,
     xlim = c(-6, 6),
     xlab = "Batch 2 Repeat",
     ylim = c(-6, 6),
     ylab = "Batch 1",
     main = "Log2 Differences in AOM-DSS vs. Control")
abline(0, 1, col = "red")

plot(dtc$log2.aomdss.cur.x ~ 
       dtc$log2.aomdss.cur.y,
     xlim = c(-6, 6),
     xlab = "Batch 2 Repeat",
     ylim = c(-6, 6),
     ylab = "Batch 1",
     main = "Log2 Differences in Curcumin vs. AOM-DSS")
abline(0, 1, col = "red")

# Absolute differences between batches----
plot(dtc$diff.aomdss.ctrl.x  ~ 
       dtc$diff.aomdss.ctrl.y,
     xlab = "Batch 2 Repeat",
     ylab = "Batch 1",
     main = "Absolute Differences in AOM-DSS vs. Control")
abline(0, 1, col = "red")

plot(dtc$diff.aomdss.cur.x ~ 
       dtc$diff.aomdss.cur.y,
     xlab = "Batch 2 Repeat",
     ylab = "Batch 1",
     main = "Absolute Differences in Curcumin vs. AOM-DSS")
abline(0, 1, col = "red")

# Plot Sample 1 vs. Sample 2 from Batch 1
x <- -log(dtc$Control_1.x)
y <- -log(dtc$Control_2.x)
  
plot(y ~ x,
     xlab = "-log2(Control) Batch 1",
     ylab = "-log2(Control) Batch 2 Rerun",
     main = "RED:log2(AOM-DSS vs. Control) <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.x) >= 1],
       col = "red")

x <- -log(dtc$AOM_DSS_1.x)
y <- -log(dtc$AOM_DSS_2.x)
plot(y ~ x,
     xlab = "-log2(AOM_DSS) Batch 1",
     ylab = "-log2(AOM_DSS) Batch 2 Rerun",
     main = "RED:log2(AOM-DSS vs. Control) <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.ctrl.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.ctrl.x) >= 1],
       col = "red")

x <- -log(dtc$Cur_1.x)
y <- -log(dtc$Cur_2.x)
plot(y ~ x,
     xlab = "-log2(Curcumin) Batch 1",
     ylab = "-log2(Curcumin) Batch 2 Rerun",
     main = "RED:log2(Curcumin vs. AOM-DSS) <= -1 or >= 1")
points(y[abs(dtc$log2.aomdss.cur.x) >= 1] ~ 
         x[abs(dtc$log2.aomdss.cur.x) >= 1],
       col = "red")

# MA Plot: X+Y vs. X-Y
M <- log2(dtc$Control_1.x) - log2(dtc$Control_2.x)
A <- (log2(dtc$Control_1.x) + log2(dtc$Control_2.x))/2
plot(M ~ A,
     main = "Batch 1 Controls")

M <- log2(dtc$AOM_DSS_1.x) - log2(dtc$AOM_DSS_2.x)
A <- (log2(dtc$AOM_DSS_1.x) + log2(dtc$AOM_DSS_2.x))/2
plot(M ~ A, 
     main = "Batch 1 AOM_DSS")

M <- log2(dtc$Cur_1.x) - log2(dtc$Cur_2.x)
A <- (log2(dtc$Cur_1.x) + log2(dtc$Cur_2.x))/2
plot(M ~ A, 
     main = "Batch 1 Curcumin")

M <- log2(dtc$aomdss_x) - log2(dtc$ctrl_x)
A <- (log2(dtc$aomdss_x) + log2(dtc$ctrl_x))/2
plot(M ~ A,
     main = "Batch 1 AOM-DSS vs. Control Averages")

M <- log2(dtc$cur_x) - log2(dtc$aomdss_x)
A <- (log2(dtc$cur_x) + log2(dtc$aomdss_x))/2
plot(M ~ A,
     main = "Batch 1 Curcumin vs. AOM-DSS Averages")

M <- log2(dtc$cur_x) - log2(dtc$ctrl_x)
A <- (log2(dtc$cur_x) + log2(dtc$ctrl_x))/2
plot(M ~ A,
     main = "Batch 1 Curcumin vs. Control Averages")

# Absolute vs. Fold-change differences----
plot(dtc$diff.aomdss.ctrl.x  ~ 
       dtc$log2.aomdss.ctrl.x,
     xlab = "log2(AOM_DSS/Control)",
     ylab = "AOM_DSS - Control",
     main = "Batch 1")

plot(dtc$diff.aomdss.ctrl.y  ~ 
       dtc$log2.aomdss.ctrl.y,
     xlab = "log2(AOM_DSS/Control)",
     ylab = "AOM_DSS - Control",
     main = "Batch 2 Rerun")


# CONTINUE HERE, DS 11/07/2017!
# DO PCA!
#   
# # Plot averages of same treatment in different batches----
# plot(dtc$ctrl_x ~ dtc$ctrl_y)
# plot(dtc$aomdss_x ~ dtc$aomdss_y)
# plot(dtc$cur_x ~ dtc$cur_y)