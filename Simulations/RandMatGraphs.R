# ____________________________________________________________________________________________
# Graphs ----
# ____________________________________________________________________________________________

# ____________________________________________________________________________________________
### NFD Maps (Fig. 2) ----
# ____________________________________________________________________________________________
ggplot(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == 120),
       aes(y = FD,
           x = ND))+
  geom_point(aes(shape = condition))+
  scale_shape_manual(values = c(1,16))+
  geom_abline(colour = "red") +
  stat_ellipse(geom = "polygon",
               fill = 4,
               colour = "black",
               alpha = 0.05,
               level = 0.95,
               linetype = 8)+
  stat_summary(aes(x = mean_ND, y = mean_FD), geom = "point", fun = "mean", colour = "black", shape = 4, size = 4)+
  ggtitle("(a)")+
  labs(y = "Fitness differences",
       x = " ")+
  xlim(0.746, 0.783)+
  ylim(0.746, 0.783)+
  theme_bw()+
  theme(legend.position = c(0.2,0.85),
        legend.background = element_rect(colour = "grey"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "57"),
       aes(y = FD,
           x = ND))+
  geom_point(aes(shape = condition))+
  scale_shape_manual(values = c(2,17))+
  geom_abline(colour = "red") +
  stat_ellipse(geom = "polygon",
               fill = 4,
               colour = "black",
               alpha = 0.05,
               level = 0.95,
               linetype = 8)+
  ylim(0.775,0.815)+
  xlim(0.775,0.815)+
  stat_summary(aes(x = mean_ND, y = mean_FD), geom = "point", fun = "mean", colour = "black", shape = 4, size = 4)+
  ggtitle("(b)")+
  labs(y = " ",
       x = " ")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "grey"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "468"),
       aes(y = FD,
           x = ND))+
  geom_point(aes(shape = condition))+
  scale_shape_manual(values = c(0,15))+
  geom_abline(colour = "red") +
  stat_ellipse(geom = "polygon",
               fill = 4,
               colour = "black",
               alpha = 0.05,
               level = 0.95,
               linetype = 8)+
  ylim(0.789,0.846)+
  xlim(0.789,0.846)+
  stat_summary(aes(x = mean_ND, y = mean_FD), geom = "point", fun = "mean", colour = "black", shape = 4, size = 4)+
  ggtitle("(c)")+
  labs(y = "Fitness differences",
       x = "Niche differences")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "grey"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = FD,
           x = ND,
           colour = lambda))+
  geom_point(alpha = 0.5,
             size = 1)+
  scale_colour_viridis_c()+
  geom_point(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "120"),
             aes(y = FD,
                 x = ND),
             shape = 1,
             size = 1,
             colour = "red")+
  stat_ellipse(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "120"),
               geom = "polygon",
               fill = 4,
               colour = "white",
               alpha = 0.05,
               level = 0.95,
               linetype = 8)+
  geom_point(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "468"),
             aes(y = FD,
                 x = ND),
             shape = 0,
             size = 1,
             colour = "red")+
  stat_ellipse(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "468"),
               geom = "polygon",
               fill = 4,
               colour = "white",
               alpha = 0.05,
               level = 0.95,
               linetype = 8)+
  geom_point(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "57"),
             aes(y = FD,
                 x = ND),
             shape = 2,
             size = 1,
             colour = "red")+
  stat_ellipse(data = subset(NFDs, condition %in% c("resident", "excluded") & matrix == "57"),
               geom = "polygon",
               fill = 4,
               colour = "white",
               alpha = 0.05,
               level = 0.95,
               linetype = 8)+
  geom_abline(colour = "red") +
  ggtitle("(d)")+
  labs(y = " ",
       x = "Niche differences",
       colour = "λ")+
  theme_bw()+
  theme(legend.position = c(0.85,0.3),
        legend.background = element_rect(colour = "grey"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# ____________________________________________________________________________________________
### Analytical vs Numerical relative differences ----
# ____________________________________________________________________________________________
ggplot(data = NFDs,
       aes(y = (c_ij - theo_c_ij)/c_ij,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(a)")+
  labs(x = " ",
       y = "Δ(A-N)/N")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (var_c_ij - theo_var_c_ij)/var_c_ij,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  # ylim(0,2)+
  ggtitle("(b)")+
  labs(x = " ",
       y = " ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_r_i - r_i)/r_i,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ylim(0,1.5)+
  ggtitle("(c)")+
  labs(x = " ",
       y = "Δ(A-N)/N")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_var_r_i - var_r_i)/var_r_i,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(d)")+
  labs(x = " ",
       y = " ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_eta_i - eta_i)/eta_i,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(e)")+
  labs(x = "Dynamics axis",
       y = "Δ(A-N)/N")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_var_eta_i - var_eta_i)/var_eta_i,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(f)")+
  labs(x = "Dynamics axis",
       y = " ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_ND - mean_ND)/mean_ND,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(a) Niche differences")+
  labs(x = " ",
       y = "Δ(A-N)/N")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_var_ND - var_ND)/var_ND,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(b) var(ND)")+
  labs(x = " ",
       y = " ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_FD - mean_FD)/mean_FD,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(c) Fitness differences")+
  labs(x = " ",
       y = "Δ(A-N)/N")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_var_FD - var_FD)/var_FD,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(d) var(FD)")+
  labs(x = " ",
       y = " ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_cov_NFD - cov_NFD)/cov_NFD,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(e) cov(ND,FD)")+
  # ylim(-0.5,3)+
  labs(x = "Dynamics axis",
       y = "Δ(A-N)/N")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(y = (theo_cor_NFD - cor_NFD)/cor_NFD,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2))))+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_hline(yintercept = 0, colour = "red")+
  stat_summary(geom = "point", fun = "mean", colour = "black", shape = 1)+
  ggtitle("(f) cor(ND,FD)")+
  # ylim(-0.5,1)+
  labs(x = "Dynamics axis",
       y = " ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# ____________________________________________________________________________________________
### Henze-Zirkler test ----
# ____________________________________________________________________________________________
ggplot(data = NFDs,
       aes(y = as.numeric(pvalue),
           x = lambda))+
  geom_point()+
  geom_hline(yintercept = 0.05, colour = "red")+
  labs(x = "Mean interaction strength (λ)",
       y = "p-value")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

p_vals = NFDs |> 
  distinct(matrix = matrix, unadj.p_vals = as.numeric(sub("^<","",pvalue))) |> 
  mutate(adj_pvals = p.adjust(unadj.p_vals, method = "BH")) |> 
  sort_by(~unadj.p_vals)

ggplot(p_vals,
       aes(x = unadj.p_vals,
           y = adj_pvals)) +
  geom_abline(slope = 1, colour = "red", linetype = "dashed") +
  geom_point() +
  ylim(0,1)+
  xlim(0,1)+
  labs(x = "Unadjusted p-value",
       y = "BH-adjusted p-value") +
  theme_bw()

# ____________________________________________________________________________________________
### Theoretical and Numerical approximations (Fig. 3) ----
# ____________________________________________________________________________________________
ggplot(data = NFDs,
       aes(x = lambda))+
  geom_point(aes(y = theo_ND), shape = 16, colour = "lightblue", size = 1)+
  geom_point(aes(y = theo_FD), shape = 15, colour = "pink", size = 1)+
  geom_smooth(aes(y = theo_ND), method = "lm", se = FALSE, colour = "blue")+
  geom_smooth(aes(y = theo_FD), method = "lm", se = FALSE, colour = "red")+
  ggtitle("(a)",
          subtitle = "Theoretical")+
  labs(y = "<ND>, <FD>",
       x = "Mean interaction strength (λ)")+
  ylim(0.73,0.84)+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(x = lambda))+
  geom_point(aes(y = mean_ND), shape = 1, colour = "lightblue", size = 1)+
  geom_point(aes(y = mean_FD), shape = 0, colour = "pink", size = 1)+
  geom_smooth(aes(y = mean_ND), method = "lm", se = FALSE, colour = "blue")+
  geom_smooth(aes(y = mean_FD), method = "lm", se = FALSE, colour = "red")+
  ggtitle("",
          subtitle = "Numerical")+
  labs(y = " ",
       x = "Mean interaction strength (λ)")+
  ylim(0.73,0.84)+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(x = rho))+
  geom_point(aes(y = theo_var_ND/((sigma**2/n)*(1-lambda)**2)), shape = 16, colour = "lightblue", size = 1)+
  geom_point(aes(y = theo_var_FD/((sigma**2/n)*(1-lambda)**2)/(1+((1-lambda)/(n*lambda)))**2), shape = 15, colour = "pink", size = 1)+
  geom_smooth(aes(y = theo_var_ND/((sigma**2/n)*(1-lambda)**2)), method = "lm", se = FALSE, colour = "blue")+
  geom_smooth(aes(y = theo_var_FD/((sigma**2/n)*(1-lambda)**2)/(1+((1-lambda)/(n*lambda)))**2), method = "lm", se = FALSE, colour = "red")+
  ylim(0,2.5)+
  ggtitle("(b)",
          subtitle = "Theoretical")+
  labs(y = "var(ND), var(FD)",
       x = "Species correlation (ρ)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))


ggplot(data = NFDs,
       aes(x = rho))+
  geom_point(aes(y = var_ND/(sigma**2/n*(1-lambda)**2)), shape = 1, colour = "lightblue", size = 1)+
  geom_point(aes(y = var_FD/(sigma**2/n*(1-lambda)**2)/(1+((1-lambda)/(n*lambda)))**2), shape = 0, colour = "pink", size = 1)+
  geom_smooth(aes(y = var_ND/(sigma**2/n*(1-lambda)**2)), method = "lm", se = FALSE, colour = "blue")+
  geom_smooth(aes(y = var_FD/(sigma**2/n*(1-lambda)**2)/(1+((1-lambda)/(n*lambda)))**2), method = "lm", se = FALSE, colour = "red")+
  ylim(0,2.5)+
  ggtitle("",
          subtitle = "Numerical")+
  labs(y = " ",
       x = "Species correlation (ρ)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(x = lambda))+
  geom_point(aes(y = cov_NFD), shape = 1, colour = "plum2", size = 1)+
  geom_point(aes(y = theo_cov_NFD), shape = 16, colour = "purple", size = 1)+
  geom_smooth(aes(y = cov_NFD), method = "lm", se = FALSE, colour = "plum2")+
  geom_smooth(aes(y = theo_cov_NFD), method = "lm", se = FALSE, colour = "purple")+
  ggtitle("(c)",
          subtitle = "Theoretical and Numerical")+
  labs(y = "cov(ND,FD)",
       x = "Mean interaction strength (λ)")+
  ylim(0, 1.5e-04)+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(x = lambda))+
  geom_point(aes(y = cor_NFD), shape = 1, colour = "plum2", size = 1)+
  geom_point(aes(y = theo_cor_NFD), shape = 16, colour = "purple", size = 1)+
  geom_smooth(aes(y = cor_NFD), method = "lm", se = FALSE, colour = "plum2")+
  geom_smooth(aes(y = theo_cor_NFD), method = "lm", se = FALSE, colour = "purple")+
  ylim(0,1)+
  ggtitle("(d)",
          subtitle = "Theoretical and Numerical")+
  labs(y = "cor(ND,FD)",
       x = "Mean interaction strength (λ)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(x = rho))+
  geom_point(aes(y = cor_NFD), shape = 1, colour = "plum2", size = 1)+
  geom_point(aes(y = theo_cor_NFD), shape = 16, colour = "purple", size = 1)+
  geom_smooth(aes(y = cor_NFD), method = "lm", se = FALSE, colour = "plum2")+
  geom_smooth(aes(y = theo_cor_NFD), method = "lm", se = FALSE, colour = "purple")+
  ylim(0,1)+
  ggtitle("(d)",
          subtitle = "Theoretical and Numerical")+
  labs(y = " ",
       x = "Species correlation (ρ)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(data = NFDs,
       aes(x = sigma))+
  geom_point(aes(y = cor_NFD), shape = 1, colour = "plum2", size = 1)+
  geom_point(aes(y = theo_cor_NFD), shape = 16, colour = "purple", size = 1)+
  geom_smooth(aes(y = cor_NFD), method = "lm", se = FALSE, colour = "plum2")+
  geom_smooth(aes(y = theo_cor_NFD), method = "lm", se = FALSE, colour = "purple")+
  ylim(0,1)+
  ggtitle("(d)",
          subtitle = "Theoretical and Numerical")+
  labs(y = " ",
       x = "Interaction standard deviation (σ)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# ____________________________________________________________________________________________
### Proportion of excluded (Fig. 4) ----
# ____________________________________________________________________________________________
ggplot(subset(NFDs, condition %in% c("excluded")),
       aes(y = proportion,
           x = (sqrt(n*sigma**2)-((1-lambda)/2)) / ((sqrt(2)/(1+rho))-((1-lambda)/2)),
           shape = condition,
           # colour = condition
       ))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = 1)+
  geom_vline(xintercept = 1, linetype = 2)+
  scale_shape_manual(values = c(16),
                     name = "Proportion",
                     labels = c("Survs / pool size"))+
  ggtitle("(a)")+
  ylim(0,1)+
  labs(y = "Proportion of excluded species",
       x = "Dynamic axis (DA)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(subset(NFDs, condition %in% c("excluded")),
       aes(y = proportion,
           x = rho,
           shape = condition,
           # colour = condition
       ))+
  geom_point()+
  ylim(0,1)+
  # scale_shape_manual(name = "Proportions",
  # labels = "Survs / pool size")+
  ggtitle("(b)")+
  labs(y = " ",
       x = "Species correlation (ρ)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

