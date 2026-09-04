# figs_plot.R -- draw Figure 1 (effect modification) and Figure 2 (multimodal concordance)
# from fig_estimates.csv. Log HR axis, reference line at 1.0, grayscale-legible (shape not
# color), estimates labeled. Saves 300-dpi PNG + vector PDF.
library(ggplot2)
est <- read.csv("fig_estimates.csv", check.names=FALSE)
## sanitize non-ASCII labels so the graphics device never has to substitute glyphs
est$label <- gsub("≤\\s*\\$?50k", "$50k or less", est$label)   # "<=$50k" -> ASCII
est$label <- gsub("≤", "<=", est$label)                        # any other <=
est$label <- gsub("≥", ">=", est$label)                        # any >=
est$label <- iconv(est$label, to = "ASCII//TRANSLIT")               # strip remaining non-ASCII
lab_ci <- function(h,l,u) sprintf("%.2f (%.2f-%.2f)", h,l,u)
pfmt   <- function(p) ifelse(p<0.001, "< .001", paste0("= ", sub("^0(\\.)", "\\1", sprintf("%.3f", p))))

## ---------- Figure 1: effect modification ----------
mod <- est[est$panel %in% c("mod_income","mod_deprivation"), ]
base <- ifelse(mod$panel=="mod_income", "A  By individual income",
               "B  By area deprivation (ZIP3)")
mod$strip <- paste0(base, "   (interaction P ", pfmt(mod$interaction_p), ")")
mod$label <- factor(mod$label, levels=rev(unique(mod$label)))

f1 <- ggplot(mod, aes(hr, label)) +
  geom_vline(xintercept=1, linetype=2, color="grey55") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0.18) +
  geom_point(size=2.8, shape=15) +
  geom_text(aes(label=lab_ci(hr,lo,hi)), hjust=0.5, vjust=-1.1, size=3) +
  facet_wrap(~strip, ncol=1, scales="free_y") +
  scale_x_log10(breaks=c(0.8,0.85,0.9,0.95,1.0)) +
  coord_cartesian(xlim=c(0.78,1.02)) +
  labs(x="Hazard ratio per 1-SD cohesion (95% CI)", y=NULL,
       title="Cohesion and incident depression\nby socioeconomic disadvantage") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(hjust=0, face="bold"),
        plot.title=element_text(size=11, face="bold"), plot.title.position="plot",
        plot.margin=margin(8,12,8,10))
ggsave("figure1_modification.png", f1, width=7.0, height=5.0, dpi=300)
ggsave("figure1_modification.pdf", f1, width=7.0, height=5.0)

## ---------- Figure 2: multimodal concordance ----------
con <- est[est$panel=="concordance", ]
grp <- c("Cohesion (exposure)"="Exposure",
         "Self-report: happiness"="Self-reported affect",
         "Self-report: meaning"="Self-reported affect",
         "Self-report: social connection"="Self-reported affect",
         "Device: steps"="Device activity","Device: MVPA"="Device activity",
         "Device: sleep efficiency"="Device sleep","Device: sleep WASO"="Device sleep")
con$grp <- factor(grp[con$label], levels=c("Exposure","Self-reported affect","Device activity","Device sleep"))
con$label <- factor(con$label, levels=rev(con$label))

f2 <- ggplot(con, aes(hr, label, shape=grp)) +
  geom_vline(xintercept=1, linetype=2, color="grey55") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0.2) +
  geom_point(size=2.8, fill="black") +
  geom_text(aes(label=lab_ci(hr,lo,hi)), hjust=0.5, vjust=-1.1, size=2.9) +
  scale_shape_manual(values=c(15,16,17,18), name=NULL) +
  scale_x_log10(breaks=c(0.6,0.7,0.8,0.9,1.0)) +
  coord_cartesian(xlim=c(0.55,1.08)) +
  labs(x="Hazard ratio per 1-SD healthier level (95% CI)", y=NULL,
       title="Association of each measurement modality\nwith incident depression",
       subtitle="All modalities protective except device-measured sleep") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        plot.title=element_text(size=11, face="bold"), plot.subtitle=element_text(size=9.5),
        plot.title.position="plot", plot.margin=margin(8,12,8,10))
ggsave("figure2_concordance.png", f2, width=7.4, height=5.2, dpi=300)
ggsave("figure2_concordance.pdf", f2, width=7.4, height=5.2)

cat("Wrote figure1_modification.{png,pdf} and figure2_concordance.{png,pdf}\n")
