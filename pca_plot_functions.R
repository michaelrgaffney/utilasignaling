pca_loadings_plot <- function(obj, components = 1:3, sortby = 1, threshold = 0, reverse=NULL){
  require(dplyr)
  if (!all(reverse %in% components)) stop("reverse vector must be a subset of components vector")
  pca_summ <- summary(obj)
  comp_nms <- colnames(pca_summ$importance)
  varprop <- pca_summ$importance[2,components] # Get % variance
  nms <- comp_nms[components]
  varprop <- paste0(nms, ' (', round(varprop*100, 1), '%)')
  names(varprop) <- nms

  pcs <- colnames(obj$rotation)[components]
  loadings <-
    obj$rotation[,components, drop = F] %*% diag(obj$sdev[components]) %>%
    data.frame %>%
    setNames(pcs)

  if(!is.null(reverse)) loadings[reverse] <- -loadings[reverse]

  loadings <- loadings[rowSums(abs(loadings)>threshold)>0,]

  loadings %>%
    dplyr::mutate(Variable = forcats::fct_reorder(rownames(.), loadings[, sortby])) %>%
    tidyr::gather(key = PC, value = Loading, -Variable) %>%
    dplyr::mutate(PC = varprop[PC]) %>%
    ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x=0, xend=Loading, y=Variable, yend=Variable, colour = Loading), size=2) +
    ggplot2::geom_point(ggplot2::aes(x=Loading, y=Variable, colour = Loading), size=3) +
    ggplot2::scale_color_gradient2(low = viridisLite::magma(11)[8], mid = 'white', 'high' = viridisLite::magma(11)[4]) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::facet_wrap(~PC) +
    ggplot2::labs(x = '\nStandardized loadings', y = '') +
    ggplot2::theme_bw(15)
}

pca_biplot <- function(obj, components = c(1,2), data = NULL, threshold = 0, reverse=NULL, label_size = 5, geom_point=T) {

  if (class(obj) != 'prcomp') stop('obj class must be prcomp')
  if (length(components) != 2 | mode(components) != 'numeric') stop('components is not a numeric vector of length 2')
  if (max(components > ncol(obj$x)) | min(components) < 1) stop(paste('components must be between 1 and ', num_pc))
  if (!all(reverse %in% components)) stop("reverse vector must be a subset of components vector")

  pcs <- colnames(obj$x)
  if (!is.null(reverse)){
    obj$x[,reverse] <- -obj$x[,reverse]
    obj$rotation[,reverse] <- -obj$rotation[,reverse]
  }

  pcX <- pcs[components[1]]
  pcY <- pcs[components[2]]

  if (!is.null(data)){
    if (!'data.frame' %in% class(data)) stop('data must be a data frame')
    if (nrow(data) != nrow(obj$x)) stop('the number of rows of data, ', nrow(data), ', does not equal the number of rows of x, ', nrow(obj$x))
    if (pcX %in% names(data) | pcY %in% names(data)) stop(paste('data cannot contain variables', pcX, 'or', pcY))
    data[pcX] <- obj$x[,pcX]
    data[pcY] <- obj$x[,pcY]
  } else {
    data = data.frame(obj$x[,components])
  }

  pct_var <- paste0('(', round(100 * summary(obj)$importance[2,components], 1), '%)')
  pct_var <- paste(colnames(obj$x)[components], pct_var)

  plot <-
    ggplot(data, aes_string(x=pcX, y=pcY)) +
    geom_hline(aes(yintercept = 0), size=.2) +
    geom_vline(aes(xintercept = 0), size=.2) +
    labs(x = pct_var[1], y = pct_var[2])

  if (geom_point){
    plot <- plot + geom_point()
  }

  datapc <- data.frame(varnames=rownames(obj$rotation), obj$rotation[,components])
  datapc <- datapc[rowSums(abs(datapc[-1])>threshold)>0,]
  mult <- min(
    (max(data[,pcY]) - min(data[,pcY])/(max(datapc[,pcY])-min(datapc[,pcY]))),
    (max(data[,pcX]) - min(data[,pcX])/(max(datapc[,pcX])-min(datapc[,pcX])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(pcX)),
                      v2 = .7 * mult * (get(pcY))
  )
  plot +
    coord_equal() +
    ggrepel::geom_text_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = label_size, color="red", max.overlaps = Inf) +
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red") +
    theme_minimal()
}
