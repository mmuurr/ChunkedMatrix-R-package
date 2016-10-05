columnwiseCosineSimilarity.chunked <- function(l, r = l, NA_as = 0,
                                               n_chunks = NULL,
                                               max_chunk_size = NULL,
                                               target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                                               mc.cores = NULL,
                                               .VERBOSE = FALSE) {
    if(.VERBOSE) timecat(sprintf("colCosSim: computing numerator... "))
    numer <- t(l) %*% r
    if(.VERBOSE) timecat(sprintf("done.\n"))
    
    if(.VERBOSE) timecat(sprintf("colCosSim: computing col norms for l... "))
    col_norms_of_l <- colNorms.chunked(m = l, p = 2,
                                       n_chunks=n_chunks,
                                       max_chunk_size=max_chunk_size,
                                       target_n_bytes_per_chunk=target_n_bytes_per_chunk, n_bytes_per_cell=n_bytes_per_cell,
                                       mc.cores=mc.cores, .VERBOSE=.VERBOSE)
    if(.VERBOSE) timecat(sprintf("done.\n"))

    if(.VERBOSE) timecat(sprintf("colCosSim: computing col norms for r... "))
    if(identical(l, r)) {
        col_norms_of_r <- col_norms_of_l
    } else {
        col_norms_of_r <- colNorms.chunked(m = r, p = 2,
                                           n_chunks=n_chunks,
                                           max_chunk_size=max_chunk_size,
                                           target_n_bytes_per_chunk=target_n_bytes_per_chunk, n_bytes_per_cell=n_bytes_per_cell,
                                           mc.cores=mc.cores, .VERBOSE=.VERBOSE)
    }
    if(.VERBOSE) timecat(sprintf("done.\n"))

    if(.VERBOSE) timecat(sprintf("colCosSim: computing numer / denom... "))
    dummy_denom <- sparseMatrix(i = integer(0), j = integer(0), x = integer(0), dims = c(ncol(l), ncol(r)))
    ## ALERT: chunking the numer, denom, and result by cols corresponds to chunking the asRow(col_norms_of_r) vector below.
    chunk_iixs <- getChunkIixs(dummy_denom, 2,
                               n_chunks=n_chunks,
                               max_chunk_size=max_chunk_size,
                               target_n_bytes_per_chunk=target_n_bytes_per_chunk, n_bytes_per_cell=n_bytes_per_cell)
    RV <- do.call(cbind,
            mclapply(seq_along(chunk_iixs), function(chunk_iixs_ix) {
                timecat(sprintf("processing chunk %d (of %d)...\n", chunk_iixs_ix, length(chunk_iixs)))
                chunk_iix <- chunk_iixs[[chunk_iixs_ix]]
                chunk_numer <- numer[,chunk_iix]
                chunk_denom <- asCol(col_norms_of_l) %*% asRow(col_norms_of_r[chunk_iix])
                Matrix(NA2(chunk_numer / chunk_denom, NA_as))
            }, mc.cores=mc.cores))
    if(.VERBOSE) timecat(sprintf("done.\n"))

    return(RV)
}
colCosSim.chunked <- columnwiseCosineSimilarity.chunked


## UNUSED.columnwiseCosineSimilarity.chunked <- function(l, r = l,
##                                                n_chunks = NULL,
##                                                max_chunk_size = NULL,
##                                                target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
##                                                mc.cores = NULL,
##                                                .VERBOSE = FALSE) {
##     numer <- t(l) %*% r
##     col_norms_of_l <- colNorms.chunked(m = l, p = 2, n_chunks = n_chunks, max_chunk_size = max_chunk_size, target_n_bytes_per_chunk = target_n_bytes_per_chunk, n_bytes_per_cell = n_bytes_per_cell, mc.cores = mc.cores, .VERBOSE = .VERBOSE)
##     if(identical(l, r)) {
##         col_norms_of_r <- col_norms_of_l
##     } else {
##         col_norms_of_r <- colNorms.chunked(m = r, p = 2, n_chunks = n_chunks, max_chunk_size = max_chunk_size, target_n_bytes_per_chunk = target_n_bytes_per_chunk, n_bytes_per_cell = n_bytes_per_cell, mc.cores = mc.cores, .VERBOSE = .VERBOSE)
##     }
##     denom <- asCol(col_norms_of_l) %*% asRow(col_norms_of_r)
##     numer / denom
## }


marginNormalize.chunked <- function(m, margin,
                                    p = 2,
                                    n_chunks = NULL,
                                    max_chunk_size = NULL,
                                    target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                                    mc.cores = NULL,
                                    .VERBOSE = FALSE) {
    v <- do.call(marginNorms.chunked, as.list(environment(), all = TRUE)[names(formals())])
    if(margin == 1) {
        Diagonal(x = 1 / v) %*% m
    } else if(margin == 2) {
        m %*% Diagonal(x = 1 / v)
    } else { stop("margin error") }
}
rowNormalize.chunked <- function(m,
                                 p = 2,
                                 n_chunks = NULL,
                                 max_chunk_size = NULL,
                                 target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                                 mc.cores = NULL,
                                 .VERBOSE = FALSE) {
    do.call(marginNormalize.chunked, c(as.list(environment(), all = TRUE)[names(formals())], margin = 1))
}
colNormalize.chunked <- function(m,
                                 p = 2,
                                 n_chunks = NULL,
                                 max_chunk_size = NULL,
                                 target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                                 mc.cores = NULL,
                                 .VERBOSE = FALSE) {
    do.call(marginNormalize.chunked, c(as.list(environment(), all = TRUE)[names(formals())], margin = 2))
}


##--------------------------------------------------------------------------------
## exactly one of `n_chunks`, `max_chunk_size`, or `target_byte_size` must be non-NULL.
## `n_bytes_per_cell` is ignored if `target_byte_size` is NULL.
##--------------------------------------------------------------------------------
marginNorms.chunked <- function(m, margin,
                                p = 2,
                                n_chunks = NULL,
                                max_chunk_size = NULL,
                                target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                                mc.cores = NULL,
                                .VERBOSE = FALSE) {
    #require(parallel)
    stopifnot(sum(!is.null(n_chunks), !is.null(max_chunk_size), !is.null(target_n_bytes_per_chunk)) == 1)
    stopifnot(length(margin) == 1 && (margin == 1 || margin == 2))
    if(is.null(mc.cores)) mc.cores <- 1
    
    if(!is.null(n_chunks)) {
        chunk_iixs <- chunk(seq_len(dim(m)[margin]), n_chunks = n_chunks, method = "seq")
    } else if(!is.null(max_chunk_size)) {
        chunk_iixs <- chunk(seq_len(dim(m)[margin]), max_chunk_size = max_chunk_size, method = "seq")
    } else if(!is.null(target_n_bytes_per_chunk)) {
        chunk_iixs <- chunk(seq_len(dim(m)[margin]), max_chunk_size = determineMaxChunkSize(m, margin = margin, target_n_bytes_per_chunk = target_n_bytes_per_chunk, n_bytes_per_cell = n_bytes_per_cell), method = "seq")
    } else { stop("argument error") }
    
    do.call(c, parallel::mclapply(seq_along(chunk_iixs), function(chunk_iixs_ix) {
        if(.VERBOSE) timecat(sprintf("processing chunk %d (of %d)...\n", chunk_iixs_ix, length(chunk_iixs)))
        if(margin == 1) {
            apply(as.matrix(m[chunk_iixs[[chunk_iixs_ix]],]), margin, lpnorm, p)
        } else if(margin == 2) {
            apply(as.matrix(m[, chunk_iixs[[chunk_iixs_ix]] ]), margin, lpnorm, p)
        } else { stop("margin error") }
    }, mc.cores = mc.cores))
}
rowNorms.chunked <- function(m,
                             p = 2,
                             n_chunks = NULL,
                             max_chunk_size = NULL,
                             target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                             mc.cores = NULL,
                             .VERBOSE = FALSE) {
    do.call(marginNorms.chunked, c(as.list(environment(), all.names = TRUE)[names(formals())], margin = 1))
}
colNorms.chunked <- function(m,
                             p = 2,
                             n_chunks = NULL,
                             max_chunk_size = NULL,
                             target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8,
                             mc.cores = NULL,
                             .VERBOSE = FALSE) {
    do.call(marginNorms.chunked, c(as.list(environment(), all.names = TRUE)[names(formals())], margin = 2))
}


determineMaxChunkSize <- function(m, margin, target_n_bytes_per_chunk, n_bytes_per_cell = 8) {
    if(margin == 1) { ## row-blocks are chunked.
        n_chunks <- max(floor(target_n_bytes_per_chunk / (ncol(m) * n_bytes_per_cell)), 1)
    } else if(margin == 2) { ## col-blocks are chunked.
        n_chunks <- max(floor(target_n_bytes_per_chunk / (nrow(m) * n_bytes_per_cell)), 1)
    }
    return(n_chunks)
}


getChunkIixs <- function(m, margin = c(1,2),
                         n_chunks = NULL,
                         max_chunk_size = NULL,
                         target_n_bytes_per_chunk = NULL, n_bytes_per_cell = 8) {
    stopifnot(length(margin) == 1 && (margin == 1 || margin == 2))
    stopifnot(sum(!is.null(n_chunks), !is.null(max_chunk_size), !is.null(target_n_bytes_per_chunk)) == 1)
    
    if(!is.null(n_chunks)) {
        return(chunk(seq_len(dim(m)[margin]), n_chunks = n_chunks, method = "seq"))
    } else if(!is.null(max_chunk_size)) {
        return(chunk(seq_len(dim(m)[margin]), max_chunk_size = max_chunk_size, method = "seq"))
    } else if(!is.null(target_n_bytes_per_chunk)) {
        return(chunk(seq_len(dim(m)[margin]),
                     max_chunk_size = determineMaxChunkSize(m, margin = margin, target_n_bytes_per_chunk = target_n_bytes_per_chunk, n_bytes_per_cell = n_bytes_per_cell),
                     method = "seq"))
    } else {
        stop("argument error")
    }
}

