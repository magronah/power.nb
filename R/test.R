# library(foreach)
# library(latex2exp)
# library(metR)
# library(scam)
# library(viridis)
# ########################################################################
# path = "/Users/michaelagronah/Documents/"
# fit_3d_0 = readRDS(paste0(path,"fit_3d_0.RDS"))
# fit_3d_1 = readRDS(paste0(path,"fit_3d_1.RDS"))
# fit_3d_2 = readRDS(paste0(path,"fit_3d_2.RDS"))
# ########################################################################
# dispersion_param  =  dispersionFit$param
# nsim = 4
# nsample_vec = seq(10, 100, 10)
# countdata_sims_list  =  list()
# for(j in 1:length(nsample_vec)){
#   countdata_sims_list[[j]]  =  countdata_sim_fun(logmean_param,
#                                                  logfoldchange_param,
#                                                  dispersion_param,
#                                                  nsamp_per_group = nsample_vec[j],
#                                                  ncont  = NULL,
#                                                  ntreat = NULL,
#                                                  notu = 1000,
#                                                  nsim = nsim,
#                                                  disp_scale = 0.3,
#                                                  max_lfc = 15,
#                                                  maxlfc_iter = 1000,
#                                                  seed = NULL)
# }
#
# names(countdata_sims_list) = paste0("sample_",nsample_vec)
# ########################################################################
# logfoldchange_list  =   read_data(countdata_sims_list,"logfoldchange_list")
# logmean_list        =   read_data(countdata_sims_list,"logmean_list")
# sample_size         =   rep(nsample_vec,
#                            times = sapply(lapply(logfoldchange_list, unlist),
#                                       length))
# ########################################################################
# dd = data.frame(logmean = unlist(logmean_list),
#                 abs_lfc = abs(unlist(logfoldchange_list)),
#                 sample_size   =  sample_size,
#                logsample_size = log(sample_size),
#                pval_reject  = sample(c(0, 1),
#                            size = length(sample_size), replace = TRUE))
# rownames(dd) = NULL
#
# ################################################################
# ddf = data.frame(mod0 = predict(fit_3d_1, newdata = dd,
#                                 type = "response"),
#                  #mod1 = predict(fit_3d_1, newdata = dd,
#                   #                     type = "response"),
#                  mod2 = predict(fit_3d_2, newdata = dd,
#                                 type = "response"))
#
#
#
# panel_diagline <- function(x, y, ...) {
#   points(x, y, ...)
#   abline(a = 0, b = 1, col = "red", lwd = 2)
# }
#
# pairs(ddf,
#       panel = panel_diagline,
#       pch = 19,
#       col = "blue")
# ################################################################
# grid_len =  20
# power_estimate  =   with(dd,
#                  expand.grid(logmean = seq(min(logmean),
#                                                max(logmean),
#                                                length  = grid_len),
#                              abs_lfc   =  seq(min(abs_lfc),
#                                               max(abs_lfc),
#                                               length  =  grid_len),
#                              sample_size = nsample_vec
#                              #logsample_size = log2(nsample_vec)
#                  #          sample_size  =  seq(min(sample_size),
#                  #                          max(sample_size),
#                  #                          length = grid_len),
#                  # logsample_size = seq(min(logsample_size),
#                  #                      max(logsample_size),
#                  #                      length = grid_len))
#                  ))
#
# power_estimate$logsample_size = log2(power_estimate$sample_size)
# power_estimate$power = predict(fit_3d_1, newdata = power_estimate,
#                                type = "response")
#
# ################################################################
# dd$pvalue_reject <- factor(dd$pval_reject)
# power_estimate$sample_size   <- factor(power_estimate$sample_size)
# ################################################################
# library(dplyr)
# sel_n <- c(10, 40, 60, 70, 100)
# power_sub <- power_estimate %>%
#   filter(sample_size %in% sel_n)
# ################################################################
# one_break = 0.2
# gg_2dimc_one <- ggplot(dd, aes(logmean, abs_lfc)) +
#   ggrastr::rasterise(
#     geom_point(color = "grey30", alpha = 0.35)
#   ) +
#   geom_contour(
#     data = power_sub,
#     aes(
#       x = logmean,
#       y = abs_lfc,
#       z = power,
#       group = sample_size,
#       colour = factor(sample_size)
#     ),
#     breaks = one_break,
#     linewidth = 1
#   ) +
#   metR::geom_label_contour(
#     data = power_sub,
#     aes(
#       x = logmean,
#       y = abs_lfc,
#       z = power,
#       colour = factor(sample_size)
#     ),
#     breaks = one_break
#   ) +
#   xlab(TeX("$\\log_2$(mean counts)")) +
#   ylab(TeX("|$\\log_2$(fold change)|")) +
#   scale_colour_discrete(name = "Sample size") +
#   theme_bw()
#
# gg_2dimc_one
# ################################################################
# uniroot_ss =  function(target_power,logmean, abs_lfc,model,xmin,xmax){
#
#   root <- uniroot(function(ss) {
#
#     data = data.frame(logsample_size = ss,
#                       logmean = logmean,
#                       abs_lfc = abs_lfc)
#
#     ##' predict.scam seem to have a bug and
#     ##' would not predict row data with only one row
#     pred = predict(model,
#                    type = "response",
#                    newdata = data[rep(1,2),])
#
#     pred[[1]] - target_power
#   },
#   interval = c(xmin, xmax),
#   extendInt = "yes")$root
#   2^root
# }
#
#
# uniroot_ss(target_power = 0.9,logmean = log2(7),
#            abs_lfc = log2(2.5),
#            model = fit_3d_0,
#            xmin = log2(0.1), xmax = log2(1000))
# ################################################################
# ##probabl use a nonlog scale
# ## try a simplier model
#
#
#
#
#
#
# cont_breaks  =  seq(0,1,0.2)
# gg_2dimc <- (ggplot(dd)
#              + aes(logmean, abs_lfc, group = sample_size)
#              + ggrastr::rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
#              + xlab(TeX("$\\log_2$(mean counts)"))
#              + ylab(TeX("|$\\log_2$(fold change)|"))
#              #+ scale_colour_manual(values = c("black", "red"))
#              + ggplot2::geom_contour(data = power_estimate,
#                                      aes(z=.data$power
#                                          # group=sample_size,
#                                          #colour = sample_size
#                                      ),lwd=1,
#                                      breaks = cont_breaks)
#              + metR::geom_label_contour(data = power_estimate,
#                                         aes(z= .data$power,label = sprintf("%.3f", after_stat(level)),
#                                             group=sample_size,
#                                             colour = sample_size
#                                         ),
#                                         breaks = cont_breaks)
#              + ggplot2::theme_bw()
#              #+ facet_wrap(~sample_size)
# )
# gg_2dimc
# ################################################################
#
#
#
#
# cont_breaks =  seq(0,1,0.1)
# cont =   contour_plot_fun(combined_data = dd,
#                           power_estimate = power_estimate,
#                           cont_breaks = cont_breaks,
#                           multiple_samples = FALSE)
#
# cont
#
#
#
#
#
# ################################################################
#
# names(power_estimate)[[1]] = "lmean_abund"
# names(power_estimate)
# names(dd)
#
#
#
# deseq_list = lapply(desq_est_list, function(x){
#   read_data(x, "deseq_estimate")
# })
#
# pval_est_list <- lapply(deseq_list, function(sample_list) {
#   lapply(sample_list, function(sim_df) {
#     sim_df$padj
#   })
# })
#
#
# logfoldchange_list  =   read_data(countdata_sims_list,"logfoldchange_list")
# logmean_list        =   read_data(countdata_sims_list,"logmean_list")
#
# ########################################################################
# fit_3d_0 <- scam::scam(pval_reject ~ s(logmean, abs_lfc,bs="tedmi") +
#                          s(logsample_size,bs="mpi"),
#                        data = comb, family = binomial)
#
#
# fit_3d_1 <- scam::scam(pval_reject ~ s(logmean, abs_lfc,bs="tedmi") +
#                          s(sample_size,logmean,bs="tedmi") +
#                          s(sample_size,abs_lfc,bs="tedmi"),
#                        data = comb, family = binomial)
#
# fit_3d_2 <- scam::scam(pval_reject ~ s(logmean, abs_lfc,bs="tedmi") +
#                          s(logsample_size,logmean,bs="tedmi") +
#                          s(logsample_size,abs_lfc,bs="tedmi"),
#                        data = comb, family = binomial)
# ########################################################################
# saveRDS(fit_3d_0, file = "/Users/michaelagronah/Documents/fit_3d_0.RDS")
# saveRDS(fit_3d, file = "/Users/michaelagronah/Documents/fit_3d_1.RDS")
# saveRDS(fit_3d_2, file = "/Users/michaelagronah/Documents/fit_3d_2.RDS")
# ########################################################################
