###
#   Author: 
#   Max Westphal
###
#   Date:
#   2022-08-24
###
#   Project:
#   cases_simstudy
###
#   Task: 
#   cases_simstudy_paper_figures.R: Reproduction of figures for main paper:
#   Westphal, Max, and Antonia Zapf. "Statistical Inference for Diagnostic Test Accuracy Studies
#   with Multiple Comparisons." arXiv preprint arXiv:2105.13469 (2021).
###



# Preparation ---------------------------------------------------------------------------------

## Load packages:
library(readr)
library(dplyr)
library(ggplot2)
library(splitstackshape)
library(glmnet)
library(cases)

## set input dir, needs to be equal to out.dir from
## cases_simstudy_lfc.R and cases_simstudy_roc.R:
in.dir <- file.path("E:/cases_simstudy/cases_sim_results")

## directory where figures are saved:
fig.dir <- file.path("figures")
dir.create(fig.dir)

## figure resolution:
fig.width <- fig.height <- 2500

## Define ggplot2 theme:
owntheme <- theme(legend.position = "bottom",
                  title=element_text(face="bold", size=24),
                  axis.title = element_text(face="bold", size=20),
                  legend.text = element_text(face="plain", size=20),
                  legend.title = element_text(face="bold", size=20), 
                  axis.text = element_text(size=16)) 

# Load data -----------------------------------------------------------------------------------

## LFC setting data:
data_lfc <- 
  readr::read_csv2(file.path(in.dir, "cases_SIM_LFC.csv"))%>%
  mutate(Adjustment = recode_factor(
    interaction(adjustment, pars),
    "none.list()" = "none",
    "bonferroni.list()"="Bonferroni", 
    "maxt.list()" = "maxT",
    "mbeta.list()" = "mBeta",
    "bootstrap.list(type='wild')" = "Bootstrap (wild)",
    "bootstrap.list()" = "Bootstrap (pairs)")
  ) %>%  
  filter(cortype == "equi") %>% 
  filter(rhose == 0.5) %>%
  filter(transformation == "none") 

dim(data_lfc)


## ROC (= 'Biomarker') setting data:
data_roc <- 
  readr::read_csv2(file.path(in.dir, "cases_SIM_ROC.csv"))%>%
  mutate(Adjustment = recode_factor(
    interaction(adjustment, pars),
    "none.list()" = "none",
    "bonferroni.list()"="Bonferroni", 
    "maxt.list()" = "maxT",
    "mbeta.list()" = "mBeta",
    "bootstrap.list(type='wild')" = "Bootstrap (wild)",
    "bootstrap.list()" = "Bootstrap (pairs)")
  ) %>% 
  filter(rhose == 0.5) %>%
  filter(auc == "seq(0.85, 0.90, length.out = 5)") %>% 
  filter(transformation == "none") 

dim(data_roc)



# LFC setting results (Figure 2) --------------------------------------------------------------
data_lfc %>% 
  filter(n<1000, prev==0.25, se==0.8, m==10) %>%
  group_by(Adjustment, n, prev, se, m) %>%
  summarize(nsim=n(), rr = mean(fp > 0)) %>%
  ggplot(aes(n, rr, col=Adjustment, fill=Adjustment, pch=Adjustment))+
  geom_point(size=6, alpha=0.6) +
  geom_line(size=2, alpha=0.6) +
  xlab("Total sample size") +
  ylab("FWER") +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  geom_hline(yintercept=0.025, col="red", lwd=2, lty=2, alpha=0.75) +
  owntheme + 
  theme(legend.position = "none")

## Figure 2.a:
ggsave("cases_sim_lfc_fwer.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)

data_lfc %>% 
  filter(n<1000, prev==0.25, se==0.8, m==10) %>%
  group_by(Adjustment, n, prev, se, m) %>%
  summarize(nsim=n(), rr = mean(`0.05` > 0)) %>%
  ggplot(aes(n, rr, col=Adjustment, fill=Adjustment, pch=Adjustment)) +
  geom_point(size=6, alpha=0.6) +
  geom_line(size=2, alpha=0.6) +
  xlab("Total sample size") +
  ylab("Power") + 
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  owntheme + 
  theme(legend.position = c(0.8, 0.4))

## Figure 2.b:
ggsave("cases_sim_lfc_power.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)

# ROC setting results (Figure 3) --------------------------------------------------------------
data_roc %>% 
  filter(n<1000, prev==0.25, m==10) %>%
  group_by(Adjustment, n, prev, m) %>%
  summarize(nsim=n(), rr = mean(`0` > 0)) %>%
  ggplot(aes(n, rr, col=Adjustment, fill=Adjustment, pch=Adjustment))+
  geom_point(size=6, alpha=0.6) +
  geom_line(size=2, alpha=0.6) +
  xlab("Total sample size") +
  ylab("FWER") +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  geom_hline(yintercept=0.025, col="red", lwd=2, lty=2, alpha=0.75) +
  owntheme +
  theme(legend.position = "none")

## Figure 3.a:
ggsave("cases_sim_roc_fwer.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)

data_roc %>% 
  filter(n<1000, prev==0.25, m==10) %>%
  group_by(Adjustment, n, prev, m) %>%
  summarize(nsim=n(), rr = mean(`0.1` > 0)) %>%
  ggplot(aes(n, rr, col=Adjustment, fill=Adjustment, pch=Adjustment)) +
  geom_point(size=6, alpha=0.6) +
  geom_line(size=2, alpha=0.6) +
  xlab("Total sample size") +
  ylab("Power") + 
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  owntheme +
  theme(legend.position = c(0.8, 0.4))

## Figure 3.b:
ggsave("cases_sim_roc_power.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)


# Synthetic example (Figure 1) ----------------------------------------------------------------

## this example is taken (copied 1:1) from the "package_overview" vignette, see
## vignette(topic="package_overview", package = "cases")

set.seed(1337)
data <- draw_data_roc(n=120, prev=c(0.25, 0.75), m=4,
                      delta=0.05, e=10, auc=seq(0.90, 0.95, 0.025), rho=c(0.25, 0.25))

head(data)

## comparison regions
results_comp <- data %>% evaluate(alternative ="greater",
                                  alpha=0.025,
                                  benchmark = c(0.7, 0.8),
                                  analysis = "co-primary",
                                  regu = TRUE,
                                  adj = "maxt")
visualize(results_comp)
## Figure 1.a:
ggsave("cases_example_syn_comp.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)



## confidence regions
results_conf <- data %>% evaluate(alternative = "greater",
                                  alpha = 0.025,
                                  benchmark = c(0.7, 0.8),
                                  analysis = "full",
                                  regu = TRUE,
                                  adj = "maxt")
visualize(results_conf)
## Figure 1.b:
ggsave("cases_example_syn_conf.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)




# Real-world data example (Figure 4) ----------------------------------------------------------

## this example is taken (reduced to necessary code) from the "example_wdbc" vignette, see
## vignette(topic="package_overview", package = "cases")

data<- opendata::load_data("wdbc") %>% 
  dplyr::mutate(diagnosis = as.factor(diagnosis)) 

sp0 <- 0.7
se0 <- 0.9
benchmark <- c(sp0, se0)

## Scenario A: Biomarker assessment
pr <- seq(0,1, 0.1) 

quantile(data$area_peak, pr) # 500, 600, 700, 800, 900 ---> area
quantile(data$compactness_peak, pr) # 0.10, 0.15, 0.20, 0.25, 0.30 ---> compactness (perimeter^2 / area - 1.0)
quantile(data$concavity_peak, pr) # 0.10, 0.15, 0.20, 0.25, 0.30 ---> concavity (severity of concave portions of the contour)

cc <- c(500, 600, 700, 800, 900, 
        0.10, 0.15, 0.20, 0.25, 0.30,
        0.10, 0.15, 0.20, 0.25, 0.30) 

comp_bm <- data %>% 
  dplyr::select(area_peak, compactness_peak, concavity_peak) %>% 
  cases::categorize(cc, rep(1:3, each=5)) %>% 
  cases::compare(labels = as.numeric(as.character(data$diagnosis)))

set.seed(1337)
results_bm <- cases::evaluate(comp_bm, benchmark = benchmark,
                              alternative = "greater", alpha = 0.025,
                              adj="boot", regu=1) 

cases::visualize(results_bm)
## Figure 4.a:
ggsave("cases_example_rwd_bm.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)

## Scenario B: Prediction model evaluation

set.seed(1337)
split <- stratified(data, c('diagnosis'), 1/3, bothSets = TRUE)
val <- split[[1]] %>% as.data.frame()
trn <- split[[2]] %>% as.data.frame()

mod <- glmnet(x=trn[,-1], y=trn[,1], family="binomial", alpha=0.25)

set.seed(1337)

aa <- c(0, 0.25, 0.5, 0.75, 1)

pred_pm <- sapply(aa, function(alpha){
  mod_pm <- cv.glmnet(x = as.matrix(trn[,-1]), y=trn[,1],
                      family = "binomial",
                      type.measure = "class",
                      alpha = alpha)
  message(paste0("cv.glmnet (alpha = ", alpha, "):"))
  print(mod_pm)
  message("+++++")
  predict(mod_pm, as.matrix(val[,-1]), type="response")
})
colnames(pred_pm) <- paste0("en", aa*100)

cc <- rep(seq(0.1, 0.5, 0.1), 5)
mm <- rep(1:5, each=5)

comp_pm <- pred_pm %>% 
  cases::categorize(cc, mm) %>% 
  cases::compare(labels = as.numeric(as.character(val$diagnosis)))

set.seed(1337)
results_pm <- cases::evaluate(comp_pm, benchmark = benchmark,
                              alternative = "greater", alpha = 0.025,
                              adj="boot", regu=1) 

cases::visualize(results_pm)
## Figure 4.b:
ggsave("cases_example_rwd_pm.jpg", path = fig.dir, width = fig.width, height = fig.height,
       units="px", device='jpg', dpi=300)





