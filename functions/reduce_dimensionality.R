do_pca <- function(ngram_comp) {
  dataset_col <- which(colnames(ngram_comp) == "Dataset")
  prcomp(ngram_comp[, c(1:(dataset_col-1), (dataset_col+1):ncol(ngram_comp))],
         center = TRUE)
}

plot_pca_res <- function(pca_res, choices = 1:2, labels) {
  ggbiplot(pca_res, choices = choices, groups = labels, ellipse = TRUE, var.axes = FALSE)
}

do_tsne <- function(ngram_comp, dims, perplexity = 30, max_iter = 1000) {
  dataset_col <- which(colnames(ngram_comp) == "Dataset")
  Rtsne(ngram_comp[, 1:(dataset_col-1)], dims = dims, perplexity = perplexity, verbose = TRUE, 
        max_iter = max_iter, check_duplicates = FALSE)
}

plot_2d_tsne <- function(tsne_res, labels) {
  ggplot(as.data.frame(tsne_res[["Y"]]), aes(x = V1, y = V2, color = labels)) +
    geom_point()
}

plot_3d_tsne <- function(tsne_res, labels, colors) {
  data.frame(tsne_res[["Y"]]) %>% 
    plot_ly(x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", color = labels, colors = colors, mode = "markers", size = 10)
}
