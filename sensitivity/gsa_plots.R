library(tidyverse); theme_set(theme_bw())

par_span <- c("total", "gb")[2]
res <- c("20ac", "9km2")[2]

resp_col <- c(pOcc="#005a32", pSB="#74c476", 
              meanN="#084594", meanNg0="#6baed6", 
              medN="#4a1486", medNg0="#9e9ac8", 
              sdN="#99000d", sdNg0="#fb6a4a")

cvDev.df <- read_csv(paste0("out/", res, "/", par_span, "/BRT_cvDev.csv"))
beta.df <- read_csv(paste0("out/", res, "/", par_span, "/BRT_betaDiv.csv"))
ri.df <- read_csv(paste0("out/", res, "/", par_span, "/BRT_RI.csv")) %>%
  mutate(response=as.factor(response))
ri.df$response <- lvls_reorder(ri.df$response, 
                               match(names(resp_col), levels(ri.df$response)))

ggplot(cvDev.df, aes(x=smp, y=Dev, group=td, colour=td)) + geom_line() +
  facet_wrap(~response, scale="free_y") + 
  ggtitle(paste0("Cross validation deviance. Parameter span = ", par_span))
ggplot(beta.df, aes(x=smp, y=beta, group=td, colour=td)) + geom_line() +
  facet_wrap(~response) + 
  ggtitle(paste0("Stability of relative influence. Parameter span = ", par_span))
ggplot(filter(ri.df, smp==max(ri.df$smp) & td==max(ri.df$td) &
                response %in% c("pOcc", "pSB", "meanN")), 
       aes(x=response, y=rel.inf, fill=response)) + coord_flip() + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") + 
  scale_fill_manual(values=resp_col) + 
  scale_x_discrete(limits=rev(c("pOcc", "pSB", "meanN"))) +
  geom_bar(stat="identity", colour="gray30") + facet_wrap(~param) + 
  ggtitle(paste0("Relative influence. Parameter span = ", par_span))

ggplot(filter(ri.df, smp==max(ri.df$smp) & td==max(ri.df$td) &
                response %in% c("pOcc", "pSB", "meanN")), 
       aes(x=param, y=rel.inf, fill=response)) + coord_flip() + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") + 
  scale_fill_manual(values=resp_col) + 
  geom_bar(stat="identity", colour="gray30") + facet_wrap(~response) + 
  ggtitle(paste0("Relative influence. Parameter span = ", par_span))


ggplot(filter(ri.df, smp==max(ri.df$smp) & td==max(ri.df$td) &
                response %in% c("pOcc", "pSB", "meanN")), 
       aes(x=response, y=rel.inf, fill=response)) + coord_flip() + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") + 
  scale_fill_manual(values=resp_col) + 
  scale_x_discrete(limits=rev(c("pOcc", "pSB", "meanN")), 
                   labels=rev(c("Pr(N)", "Pr(B)", "mean(N)"))) +
  scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), 
                     labels=c("0", "0.5", "1")) + 
  geom_bar(stat="identity", colour="gray30") + facet_wrap(~param) + 
  ggtitle("Preliminary sensitivity analysis") + 
  theme(legend.position="none", 
        axis.text=element_text(size=12),
        title=element_text(size=16),
        strip.text=element_text(size=12)) + 
  labs(y="Relative Influence", x="Response Metric")







