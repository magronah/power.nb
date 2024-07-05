# read_data <- function(dataset_list, extract_name){
#   extract_data_list= list()
#   for(n in 1:length(dataset_list)){
#     dataset <- dataset_list[[n]]
#     extract_data_list[[n]] <- dataset[[extract_name]]
#   }
#   names(extract_data_list) <- names(dataset_list)
#   extract_data_list
# }
# ###############################################################################
# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
#
# ###############################################################################
# custom_theme <- function(n) {
#   theme_bw(base_size = n) +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       text = element_text(size = n, family = "Roboto"),
#       axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
#       axis.text.y = element_text(family = "Roboto", size = n, color = "black")
#     )
# }
#
# ####################################
# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
#
# initial_cond_est <- function(p0,logmean, logfoldchange,
#                              np, sd_ord, ncore,
#                              minval = -5, maxval = 5,
#                              NP = 800, itermax = 1500){
#
#   cl <- makeCluster(ncore)
#   fun_names <- c("polyfun", "dnormmix", "dnormmix0", "genmixpars")
#   clusterExport(cl, c("logmean", "logfoldchange", fun_names))
#
#   opt  <- DEoptim(fn = nllfun, lower = rep(minval, length(p0)), upper = rep(maxval, length(p0)),
#                   vals = logfoldchange, logmean = logmean, np = np, sd_ord = sd_ord,
#                   control = DEoptim.control(NP = NP, itermax = itermax,
#                                             cluster = cl))
#
#   opt$optim$bestmem
#
# }
