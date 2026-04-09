# =========================================================
# Mfuzz pipeline with YAML config (CSV input)
# ---------------------------------------------------------
# 사용법 1: RStudio에서 Source 버튼
# - config.yaml 과 input.csv 가 같은 폴더에 있으면 바로 실행됨
#
# 사용법 2: Console에서
#   source("mfuzz_pipeline_csv.R")
#
# 사용법 3: 터미널에서
#   Rscript mfuzz_pipeline_csv.R config.yaml
# =========================================================


# -----------------------------
# [0] helper functions
# -----------------------------

# NULL이면 오른쪽 값을 반환하는 helper
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# 숫자형 변환 helper
safe_as_numeric <- function(x) {
  if (is.null(x)) return(NULL)
  as.numeric(x)
}

# 파일명에 안전한 문자열로 변환
file_safe <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", tolower(x))
}

# csv 저장 helper
write_csv_safe <- function(x, file) {
  write.csv(x, file = file, row.names = FALSE, quote = FALSE, na = "")
}


# -----------------------------
# [1] install / load packages
# -----------------------------

# CRAN 패키지
required_cran <- c("yaml", "viridisLite")

# Bioconductor 패키지
required_bioc <- c("Biobase", "Mfuzz")

# BiocManager가 없으면 설치
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN 패키지 설치
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor 패키지 설치
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# 라이브러리 로드
library(yaml)
library(viridisLite)
library(Biobase)
library(Mfuzz)


# -----------------------------
# [2] palette resolver
# -----------------------------
# config에서 지정한 palette 이름을 실제 색 벡터로 변환
resolve_palette <- function(name, n = 64, reverse = FALSE) {
  palette_name <- tolower(name)
  
  # default 스타일용 수동 palette
  default_pal <- c(
    "#FF8F00", "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00",
    "#DFFF00", "#C7FF00", "#AFFF00", "#97FF00", "#80FF00", "#68FF00",
    "#50FF00", "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28",
    "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7",
    "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", "#00B7FF",
    "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF",
    "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF",
    "#8000FF", "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF",
    "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078",
    "#FF0060", "#FF0048", "#FF0030", "#FF0018"
  )
  
  # fancy 스타일용 palette
  fancy_blue  <- c(255:0, rep(0, length(255:0)), rep(0, length(255:150)))
  fancy_green <- c(0:255, 255:0, rep(0, length(255:150)))
  fancy_red   <- c(0:255, rep(255, length(255:0)), 255:150)
  fancy_pal   <- rgb(r = fancy_red / 255, g = fancy_green / 255, b = fancy_blue / 255)
  
  # 이름별 palette 반환
  pal <- switch(
    palette_name,
    "default" = grDevices::colorRampPalette(default_pal)(n),
    "fancy"   = grDevices::colorRampPalette(fancy_pal)(n),
    "viridis" = viridisLite::viridis(n),
    "magma"   = viridisLite::magma(n),
    "inferno" = viridisLite::inferno(n),
    "plasma"  = viridisLite::plasma(n),
    "cividis" = viridisLite::cividis(n),
    "hot_cold" = grDevices::colorRampPalette(
      c("#3B4CC0", "#8DB0FE", "#F7F7F7", "#F4987A", "#B40426")
    )(n),
    viridisLite::viridis(n)
  )
  
  # palette 순서 뒤집기
  if (isTRUE(reverse)) {
    pal <- rev(pal)
  }
  
  pal
}


# -----------------------------
# [3] config validation
# -----------------------------
# config.yaml 필수 항목 검증
validate_config <- function(cfg) {
  if (is.null(cfg$data)) stop("config.yaml에 data 섹션이 없습니다.")
  if (is.null(cfg$preprocess)) stop("config.yaml에 preprocess 섹션이 없습니다.")
  if (is.null(cfg$clustering)) stop("config.yaml에 clustering 섹션이 없습니다.")
  if (is.null(cfg$plot)) stop("config.yaml에 plot 섹션이 없습니다.")
  if (is.null(cfg$output)) stop("config.yaml에 output 섹션이 없습니다.")
  
  if (is.null(cfg$data$input_file)) stop("data.input_file 이 필요합니다.")
  if (is.null(cfg$data$gene_id_col)) stop("data.gene_id_col 이 필요합니다.")
  if (is.null(cfg$data$sample_groups) || length(cfg$data$sample_groups) < 2) {
    stop("data.sample_groups 는 최소 2개 이상 정의해야 합니다.")
  }
  
  for (g in cfg$data$sample_groups) {
    if (is.null(g$name) || is.null(g$time) || is.null(g$columns)) {
      stop("각 sample_group에는 name, time, columns 가 필요합니다.")
    }
    if (length(g$columns) < 1) {
      stop("각 sample_group의 columns 는 최소 1개 이상이어야 합니다.")
    }
  }
  
  if (isTRUE(cfg$clustering$estimate_m) == FALSE && is.null(cfg$clustering$m)) {
    stop("estimate_m=false 이면 clustering.m 값을 지정해야 합니다.")
  }
}


# -----------------------------
# [4] read input data
# -----------------------------
# 입력 CSV/TSV 파일을 읽고 expression matrix와 gene annotation 생성
read_expression_data <- function(cfg) {
  input_file <- cfg$data$input_file
  header     <- cfg$data$header %||% TRUE
  delimiter  <- cfg$data$delimiter %||% ","
  na_strings <- unlist(cfg$data$na_strings %||% c("NA", "", "NaN"))
  
  dat <- read.table(
    file = input_file,
    header = header,
    sep = delimiter,
    quote = "\"",
    na.strings = na_strings,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  
  gene_id_col <- cfg$data$gene_id_col
  symbol_col  <- cfg$data$symbol_col
  
  # sample_groups에서 샘플 컬럼명 수집
  sample_groups <- cfg$data$sample_groups
  sample_cols <- unlist(lapply(sample_groups, function(g) unlist(g$columns)))
  
  # 필수 컬럼명 확인
  required_cols <- c(gene_id_col, sample_cols)
  if (!is.null(symbol_col)) {
    required_cols <- c(required_cols, symbol_col)
  }
  
  if (!all(required_cols %in% colnames(dat))) {
    missing_cols <- setdiff(required_cols, colnames(dat))
    stop(
      paste0(
        "입력 파일에 필요한 컬럼이 없습니다: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # gene_id / symbol 추출
  gene_id <- as.character(dat[[gene_id_col]])
  symbol  <- if (!is.null(symbol_col)) as.character(dat[[symbol_col]]) else rep(NA_character_, nrow(dat))
  
  # 중복 gene_id가 있을 수 있으므로 내부 row_id는 unique하게 생성
  row_id <- make.unique(gene_id)
  
  # expression matrix 추출
  expr_mat <- as.matrix(dat[, sample_cols, drop = FALSE])
  rownames(expr_mat) <- row_id
  storage.mode(expr_mat) <- "numeric"
  
  # gene annotation table 생성
  gene_map <- data.frame(
    row_id = row_id,
    gene_id = gene_id,
    symbol = symbol,
    stringsAsFactors = FALSE
  )
  
  list(expr_mat = expr_mat, gene_map = gene_map)
}


# -----------------------------
# [5] aggregate replicates
# -----------------------------
# 각 시간대 replicate를 mean 또는 median으로 묶음
aggregate_groups <- function(expr_mat, sample_groups, aggregate_fun = "mean") {
  agg_one <- function(mat, cols, fun_name) {
    sub <- mat[, cols, drop = FALSE]
    
    if (fun_name == "mean") {
      out <- apply(sub, 1, function(x) {
        if (all(is.na(x))) return(NA_real_)
        mean(x, na.rm = TRUE)
      })
    } else if (fun_name == "median") {
      out <- apply(sub, 1, function(x) {
        if (all(is.na(x))) return(NA_real_)
        median(x, na.rm = TRUE)
      })
    } else {
      stop("aggregate_fun 은 mean 또는 median 만 지원합니다.")
    }
    
    out
  }
  
  group_names <- vapply(sample_groups, function(g) as.character(g$name), character(1))
  time_points <- vapply(sample_groups, function(g) as.numeric(g$time), numeric(1))
  
  group_list <- lapply(sample_groups, function(g) {
    agg_one(expr_mat, cols = unlist(g$columns), fun_name = aggregate_fun)
  })
  
  group_mat <- do.call(cbind, group_list)
  colnames(group_mat) <- group_names
  rownames(group_mat) <- rownames(expr_mat)
  
  list(group_mat = group_mat, group_names = group_names, time_points = time_points)
}


# -----------------------------
# [6] attach gene annotation
# -----------------------------
# matrix/data.frame 앞에 gene_id, symbol 컬럼을 붙임
attach_gene_cols <- function(mat_or_df, gene_map, row_ids = NULL) {
  if (is.null(row_ids)) {
    row_ids <- rownames(mat_or_df)
  }
  
  out_df <- as.data.frame(mat_or_df, check.names = FALSE, stringsAsFactors = FALSE)
  out_df$gene_id <- gene_map$gene_id[match(row_ids, gene_map$row_id)]
  out_df$symbol  <- gene_map$symbol[match(row_ids, gene_map$row_id)]
  
  out_df <- out_df[, c("gene_id", "symbol", setdiff(colnames(out_df), c("gene_id", "symbol"))), drop = FALSE]
  out_df
}


# -----------------------------
# [7] build ExpressionSet
# -----------------------------
# Mfuzz용 ExpressionSet 생성
build_eset <- function(mat, group_names, time_points) {
  pheno_df <- data.frame(
    sample = group_names,
    time = time_points,
    group = group_names,
    row.names = group_names,
    stringsAsFactors = FALSE
  )
  
  pheno_annot <- new("AnnotatedDataFrame", data = pheno_df)
  
  ExpressionSet(
    assayData = as.matrix(mat),
    phenoData = pheno_annot
  )
}


# -----------------------------
# [8] custom plot helper
# -----------------------------
# membership 값을 palette index로 변환
map_membership_to_color <- function(mem, palette_vec) {
  mem <- pmax(0, pmin(1, mem))
  idx <- floor(mem * (length(palette_vec) - 1)) + 1
  palette_vec[idx]
}

# custom colorbar 그리기
plot_colorbar_custom <- function(palette_vec, horizontal = FALSE, main = "Membership", cex.main = 1) {
  n <- length(palette_vec)
  
  if (horizontal) {
    par(mar = c(4, 2, 3, 2))
    plot(
      NA, xlim = c(0, 1), ylim = c(0, 1),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
      main = main, cex.main = cex.main
    )
    xbreaks <- seq(0, 1, length.out = n + 1)
    for (i in seq_len(n)) {
      rect(xbreaks[i], 0, xbreaks[i + 1], 1, col = palette_vec[i], border = NA)
    }
    axis(1, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
  } else {
    par(mar = c(4, 2, 3, 4))
    plot(
      NA, xlim = c(0, 1), ylim = c(0, 1),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
      main = main, cex.main = cex.main
    )
    ybreaks <- seq(0, 1, length.out = n + 1)
    for (i in seq_len(n)) {
      rect(0, ybreaks[i], 1, ybreaks[i + 1], col = palette_vec[i], border = NA)
    }
    axis(4, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2), las = 1)
  }
}

# cluster 하나를 그리는 custom plot
plot_cluster_custom <- function(
    eset_std,
    cl,
    cluster_id,
    min_mem,
    palette_vec,
    cfg,
    main_title = NULL
) {
  expr_mat <- exprs(eset_std)
  cluster_index <- cl$cluster
  membership_mat <- cl$membership
  
  keep_idx <- which(cluster_index == cluster_id)
  
  plot_extra <- cfg$plot$plot_extra %||% list()
  cex_main <- plot_extra$cex.main %||% 1.0
  cex_lab  <- plot_extra$cex.lab %||% 1.0
  cex_axis <- plot_extra$cex.axis %||% 0.9
  
  # 비어 있는 cluster면 빈 plot 출력
  if (length(keep_idx) == 0) {
    plot.new()
    title(main = paste("Cluster", cluster_id, "(empty)"), cex.main = cex_main)
    return(invisible(NULL))
  }
  
  tmp <- expr_mat[keep_idx, , drop = FALSE]
  tmpmem <- membership_mat[keep_idx, cluster_id]
  
  # min_mem 이상인 유전자만 plot
  keep_mem <- tmpmem >= min_mem
  tmp <- tmp[keep_mem, , drop = FALSE]
  tmpmem <- tmpmem[keep_mem]
  
  if (nrow(tmp) == 0) {
    plot.new()
    title(main = paste("Cluster", cluster_id, "(no genes above min.mem)"), cex.main = cex_main)
    return(invisible(NULL))
  }
  
  # membership 낮은 것 -> 높은 것 순서로 그리면 core gene이 위에 보임
  ord <- order(tmpmem)
  tmp <- tmp[ord, , drop = FALSE]
  tmpmem <- tmpmem[ord]
  
  xvals <- pData(eset_std)$time
  xlabels <- sampleNames(eset_std)
  
  ylim_set <- safe_as_numeric(cfg$plot$ylim_set %||% c(0, 0))
  if (all(ylim_set == c(0, 0))) {
    ylim_use <- range(tmp, na.rm = TRUE)
  } else {
    ylim_use <- ylim_set
  }
  
  plot(
    NA,
    xlim = c(min(xvals), max(xvals)),
    ylim = ylim_use,
    xlab = cfg$plot$xlab %||% "Time",
    ylab = cfg$plot$ylab %||% "Standardized expression",
    main = main_title %||% paste("Cluster", cluster_id),
    axes = FALSE,
    bg = cfg$plot$bg %||% "white",
    col.lab = cfg$plot$col_lab %||% "black",
    col.main = cfg$plot$col_main %||% "black",
    col.sub = cfg$plot$col_sub %||% "black",
    cex.main = cex_main,
    cex.lab = cex_lab
  )
  
  axis(
    1,
    at = xvals,
    labels = xlabels,
    col = cfg$plot$ax_col %||% "black",
    col.axis = cfg$plot$col_axis %||% "black",
    cex.axis = cex_axis
  )
  axis(
    2,
    col = cfg$plot$ax_col %||% "black",
    col.axis = cfg$plot$col_axis %||% "black",
    cex.axis = cex_axis
  )
  box(col = cfg$plot$ax_col %||% "black")
  
  line_cols <- map_membership_to_color(tmpmem, palette_vec)
  
  # 각 gene profile 선 그리기
  for (i in seq_len(nrow(tmp))) {
    lines(xvals, tmp[i, ], col = line_cols[i])
  }
  
  # cluster 중심선 그리기
  if (isTRUE(cfg$plot$centre %||% TRUE)) {
    lines(
      xvals,
      cl$centers[cluster_id, ],
      col = cfg$plot$centre_col %||% "black",
      lwd = cfg$plot$centre_lwd %||% 2
    )
  }
}

# PDF 전체 plot + colorbar + cluster별 PNG 생성
render_custom_plot_set <- function(
    eset_std,
    cl,
    cfg,
    palette_name,
    min_mem,
    pdf_file,
    png_prefix
) {
  palette_vec <- resolve_palette(
    name = palette_name,
    n = cfg$plot$palette_n %||% 64,
    reverse = cfg$plot$reverse_palette %||% FALSE
  )
  
  all_cluster_ids <- seq_len(nrow(cl$centers))
  single_cluster <- cfg$plot$single_cluster
  
  cluster_ids <- if (is.null(single_cluster)) {
    all_cluster_ids
  } else {
    as.numeric(single_cluster)
  }
  
  cluster_num_for_plot <- length(cluster_ids)
  
  ncol_cfg <- cfg$plot$ncol
  ncol_plot <- if (is.null(ncol_cfg)) {
    if (cluster_num_for_plot <= 4) 2 else 3
  } else {
    as.numeric(ncol_cfg)
  }
  nrow_plot <- ceiling(cluster_num_for_plot / ncol_plot)
  
  # PDF 저장
  pdf(
    file = pdf_file,
    width = cfg$output$pdf$width %||% 12,
    height = cfg$output$pdf$height %||% 8
  )
  
  par(mfrow = c(nrow_plot, ncol_plot))
  
  for (cluster_id in cluster_ids) {
    plot_cluster_custom(
      eset_std = eset_std,
      cl = cl,
      cluster_id = cluster_id,
      min_mem = min_mem,
      palette_vec = palette_vec,
      cfg = cfg,
      main_title = paste0("Cluster ", cluster_id)
    )
  }
  
  dev.off()
  
  # colorbar 저장
  if (isTRUE(cfg$plot$colorbar$enabled %||% TRUE)) {
    colorbar_file <- sub("\\.pdf$", paste0("__", file_safe(palette_name), "_colorbar.pdf"), pdf_file)
    
    pdf(file = colorbar_file, width = 2.5, height = 5.5)
    plot_colorbar_custom(
      palette_vec = palette_vec,
      horizontal = cfg$plot$colorbar$horizontal %||% FALSE,
      main = paste(cfg$plot$colorbar$main_prefix %||% "Membership", "-", palette_name),
      cex.main = cfg$plot$colorbar$cex.main %||% 1.0
    )
    dev.off()
  }
  
  # cluster별 PNG 저장
  if (isTRUE(cfg$output$png$save_cluster_png %||% TRUE)) {
    for (cluster_id in cluster_ids) {
      png_file <- paste0(
        png_prefix,
        "__",
        file_safe(palette_name),
        "__cluster_",
        cluster_id,
        ".png"
      )
      
      png(
        filename = png_file,
        width = cfg$output$png$width %||% 1800,
        height = cfg$output$png$height %||% 1400,
        res = cfg$output$png$res %||% 200,
        bg = cfg$output$png$bg %||% "white"
      )
      
      plot_cluster_custom(
        eset_std = eset_std,
        cl = cl,
        cluster_id = cluster_id,
        min_mem = min_mem,
        palette_vec = palette_vec,
        cfg = cfg,
        main_title = paste0("Cluster ", cluster_id, " - ", palette_name)
      )
      
      dev.off()
    }
  }
}


# -----------------------------
# [9] main pipeline
# -----------------------------
run_mfuzz_pipeline <- function(config_file = "config.yaml") {
  # config 읽기
  cfg <- yaml::read_yaml(config_file)
  validate_config(cfg)
  
  # output 폴더 생성
  outdir <- cfg$output$outdir %||% "mfuzz_out"
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # 입력 데이터 읽기
  dat_obj <- read_expression_data(cfg)
  expr_mat <- dat_obj$expr_mat
  gene_map <- dat_obj$gene_map
  
  # 음수 값 체크
  if (any(expr_mat < 0, na.rm = TRUE)) {
    stop("입력값에 음수가 있습니다. 파일을 확인하세요.")
  }
  
  # 전부 NA인 gene 제거
  keep_not_all_na <- rowSums(!is.na(expr_mat)) > 0
  expr_mat <- expr_mat[keep_not_all_na, , drop = FALSE]
  gene_map <- gene_map[match(rownames(expr_mat), gene_map$row_id), , drop = FALSE]
  
  message("Initial gene count: ", nrow(expr_mat))
  
  # replicate를 시간대별로 묶기
  agg <- aggregate_groups(
    expr_mat = expr_mat,
    sample_groups = cfg$data$sample_groups,
    aggregate_fun = cfg$preprocess$aggregate_fun %||% "mean"
  )
  
  group_mat <- agg$group_mat
  group_names <- agg$group_names
  time_points <- agg$time_points
  
  # 필요시 log2 변환
  if (isTRUE(cfg$preprocess$apply_log2 %||% FALSE)) {
    pseudo_count <- cfg$preprocess$pseudo_count %||% 1
    group_mat <- log2(group_mat + pseudo_count)
  }
  
  # 중간 결과 저장
  if (isTRUE(cfg$output$save_intermediate_csv %||% TRUE)) {
    write_csv_safe(
      attach_gene_cols(group_mat, gene_map),
      file.path(outdir, "01_group_matrix_used_for_mfuzz.csv")
    )
  }
  
  # ExpressionSet 생성
  eset <- build_eset(group_mat, group_names, time_points)
  
  # 결측 비율이 너무 높은 gene 제거
  if (isTRUE(cfg$preprocess$filter_NA$enabled %||% TRUE)) {
    eset <- filter.NA(
      eset,
      thres = as.numeric(cfg$preprocess$filter_NA$thres %||% 0.34)
    )
  }
  
  # 남은 결측값 보정
  if (isTRUE(cfg$preprocess$fill_NA$enabled %||% TRUE)) {
    if (any(is.na(exprs(eset)))) {
      eset <- fill.NA(
        eset,
        mode = as.character(cfg$preprocess$fill_NA$mode %||% "mean"),
        k = as.numeric(cfg$preprocess$fill_NA$k %||% 10)
      )
    }
  }
  
  # NA 처리 후 저장
  if (isTRUE(cfg$output$save_intermediate_csv %||% TRUE)) {
    write_csv_safe(
      attach_gene_cols(exprs(eset), gene_map),
      file.path(outdir, "02_after_na_handling.csv")
    )
  }
  
  # 변동성이 낮은 gene 제거
  if (isTRUE(cfg$preprocess$filter_std$enabled %||% TRUE)) {
    min_std <- as.numeric(cfg$preprocess$filter_std$min_std %||% 0.5)
    visu <- isTRUE(cfg$preprocess$filter_std$visu %||% FALSE)
    
    if (visu) {
      pdf(
        file.path(outdir, "03_filter_std_plot.pdf"),
        width = cfg$output$pdf$width %||% 12,
        height = cfg$output$pdf$height %||% 8
      )
      eset <- filter.std(eset, min.std = min_std, visu = TRUE)
      dev.off()
    } else {
      eset <- filter.std(eset, min.std = min_std, visu = FALSE)
    }
  }
  
  # 너무 적게 남으면 중단
  if (nrow(exprs(eset)) < 10) {
    stop("필터 후 남은 gene 수가 너무 적습니다. config 값을 완화하세요.")
  }
  
  # SD filter 후 저장
  if (isTRUE(cfg$output$save_intermediate_csv %||% TRUE)) {
    write_csv_safe(
      attach_gene_cols(exprs(eset), gene_map),
      file.path(outdir, "04_after_sd_filter.csv")
    )
  }
  
  # gene-wise standardisation
  eset_std <- if (isTRUE(cfg$preprocess$standardise$enabled %||% TRUE)) {
    standardise(eset)
  } else {
    eset
  }
  
  # standardised matrix 저장
  write_csv_safe(
    attach_gene_cols(exprs(eset_std), gene_map),
    file.path(outdir, "05_zscore_or_standardised_matrix.csv")
  )
  
  # m 값 추정 또는 수동 사용
  m_val <- if (isTRUE(cfg$clustering$estimate_m %||% TRUE)) {
    mestimate(eset_std)
  } else {
    as.numeric(cfg$clustering$m)
  }
  
  write_csv_safe(
    data.frame(parameter = "m", value = m_val),
    file.path(outdir, "06_mestimate.csv")
  )
  
  # cselection 실행
  if (isTRUE(cfg$parameter_search$run_cselection %||% FALSE)) {
    c_cfg <- cfg$parameter_search$cselection
    crange <- as.numeric(unlist(c_cfg$crange %||% c(3, 4, 5, 6, 7, 8)))
    repeats <- as.numeric(c_cfg$repeats %||% 5)
    visu <- isTRUE(c_cfg$visu %||% TRUE)
    
    if (visu) {
      pdf(
        file.path(outdir, "06A_cselection_plot.pdf"),
        width = cfg$output$pdf$width %||% 12,
        height = cfg$output$pdf$height %||% 8
      )
      csel <- cselection(eset_std, m = m_val, crange = crange, repeats = repeats, visu = TRUE)
      dev.off()
    } else {
      csel <- cselection(eset_std, m = m_val, crange = crange, repeats = repeats, visu = FALSE)
    }
    
    csel_df <- data.frame(c = rownames(csel), csel, check.names = FALSE)
    write_csv_safe(csel_df, file.path(outdir, "06A_cselection_matrix.csv"))
  }
  
  # Dmin 실행
  if (isTRUE(cfg$parameter_search$run_Dmin %||% FALSE)) {
    d_cfg <- cfg$parameter_search$Dmin
    crange <- as.numeric(unlist(d_cfg$crange %||% c(3, 4, 5, 6, 7, 8)))
    repeats <- as.numeric(d_cfg$repeats %||% 3)
    visu <- isTRUE(d_cfg$visu %||% TRUE)
    
    if (visu) {
      pdf(
        file.path(outdir, "06B_Dmin_plot.pdf"),
        width = cfg$output$pdf$width %||% 12,
        height = cfg$output$pdf$height %||% 8
      )
      dmin <- Dmin(eset_std, m = m_val, crange = crange, repeats = repeats, visu = TRUE)
      dev.off()
    } else {
      dmin <- Dmin(eset_std, m = m_val, crange = crange, repeats = repeats, visu = FALSE)
    }
    
    dmin_df <- data.frame(c = crange, Dmin = dmin)
    write_csv_safe(dmin_df, file.path(outdir, "06B_Dmin_values.csv"))
  }
  
  # 최종 Mfuzz clustering 수행
  set.seed(as.integer(cfg$clustering$random_seed %||% 1234))
  
  cl <- mfuzz(
    eset = eset_std,
    c = as.numeric(cfg$clustering$cluster_number %||% 4),
    m = m_val
  )
  
  # membership matrix와 대표 cluster 계산
  membership_mat <- cl$membership
  best_cluster <- apply(membership_mat, 1, which.max)
  max_membership <- apply(membership_mat, 1, max)
  
  # gene별 대표 cluster 저장
  cluster_out <- data.frame(
    gene_id = gene_map$gene_id[match(rownames(membership_mat), gene_map$row_id)],
    symbol = gene_map$symbol[match(rownames(membership_mat), gene_map$row_id)],
    cluster = best_cluster,
    max_membership = max_membership,
    stringsAsFactors = FALSE
  )
  write_csv_safe(cluster_out, file.path(outdir, "07_gene_cluster_assignment.csv"))
  
  # alpha-core gene 추출
  alpha_cut <- as.numeric(cfg$plot$min_mem_core %||% 0.7)
  core_list <- acore(eset_std, cl = cl, min.acore = alpha_cut)
  
  total_alpha_core_list <- list()
  
  for (i in seq_along(core_list)) {
    core_i <- core_list[[i]]
    
    if (!is.null(core_i) && nrow(core_i) > 0) {
      core_df <- data.frame(
        row_id = rownames(core_i),
        core_i,
        stringsAsFactors = FALSE
      )
      
      core_df$gene_id <- gene_map$gene_id[match(core_df$row_id, gene_map$row_id)]
      core_df$symbol  <- gene_map$symbol[match(core_df$row_id, gene_map$row_id)]
      core_df$cluster_id <- i
      
      # cluster별 alpha-core 파일 저장
      core_df_out <- core_df[, c(
        "gene_id", "symbol",
        setdiff(colnames(core_df), c("gene_id", "symbol"))
      ), drop = FALSE]
      
      write_csv_safe(
        core_df_out,
        file.path(outdir, paste0("08_alpha_core_cluster_", i, ".csv"))
      )
      
      # membership 컬럼 자동 탐지
      mem_candidates <- colnames(core_df)[grepl("mem", colnames(core_df), ignore.case = TRUE)]
      if (length(mem_candidates) == 0) {
        stop("acore 결과에서 membership 컬럼명을 찾지 못했습니다.")
      }
      
      mem_col <- mem_candidates[1]
      
      # total alpha-core용 축적
      total_alpha_core_list[[length(total_alpha_core_list) + 1]] <- data.frame(
        gene_id = core_df$gene_id,
        symbol = core_df$symbol,
        MEM.SHIP = core_df[[mem_col]],
        cluster_id = i,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # alpha-core 전체 합친 파일 저장
  if (length(total_alpha_core_list) > 0) {
    alpha_core_total <- do.call(rbind, total_alpha_core_list)
  } else {
    alpha_core_total <- data.frame(
      gene_id = character(0),
      symbol = character(0),
      MEM.SHIP = numeric(0),
      cluster_id = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  write_csv_safe(
    alpha_core_total,
    file.path(outdir, "08_alpha_core_cluster_total.csv")
  )
  
  # plot palette 목록 결정
  palette_names <- if (isTRUE(cfg$plot$render_multiple_palettes %||% FALSE)) {
    unlist(cfg$plot$palette_names_to_render %||% "viridis")
  } else {
    cfg$plot$active_palette %||% "viridis"
  }
  
  palette_names <- unique(as.character(palette_names))
  
  # palette별 PDF / PNG 생성
  for (pal_name in palette_names) {
    pal_slug <- file_safe(pal_name)
    
    render_custom_plot_set(
      eset_std = eset_std,
      cl = cl,
      cfg = cfg,
      palette_name = pal_name,
      min_mem = as.numeric(cfg$plot$min_mem_all %||% 0.0),
      pdf_file = file.path(outdir, paste0("09_mfuzz_plot_all_genes__", pal_slug, ".pdf")),
      png_prefix = file.path(outdir, "09_mfuzz_plot_all_genes")
    )
    
    render_custom_plot_set(
      eset_std = eset_std,
      cl = cl,
      cfg = cfg,
      palette_name = pal_name,
      min_mem = as.numeric(cfg$plot$min_mem_core %||% 0.7),
      pdf_file = file.path(outdir, paste0("10_mfuzz_plot_core_genes__", pal_slug, ".pdf")),
      png_prefix = file.path(outdir, "10_mfuzz_plot_core_genes")
    )
  }
  
  # membership matrix 저장
  write_csv_safe(
    attach_gene_cols(membership_mat, gene_map),
    file.path(outdir, "11_membership_matrix.csv")
  )
  
  # cluster center 저장
  centers_out <- data.frame(
    cluster = seq_len(nrow(cl$centers)),
    cl$centers,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_csv_safe(centers_out, file.path(outdir, "12_cluster_centers.csv"))
  
  # cluster별 gene 수 요약 저장
  cluster_size <- as.data.frame(table(cluster_out$cluster))
  colnames(cluster_size) <- c("cluster", "gene_count")
  write_csv_safe(cluster_size, file.path(outdir, "13_cluster_size_summary.csv"))
  
  # 실행 요약 저장
  summary_out <- data.frame(
    metric = c(
      "genes_initial",
      "genes_after_final_filter",
      "cluster_number",
      "estimated_m"
    ),
    value = c(
      nrow(expr_mat),
      nrow(exprs(eset_std)),
      as.numeric(cfg$clustering$cluster_number %||% 4),
      m_val
    ),
    stringsAsFactors = FALSE
  )
  write_csv_safe(summary_out, file.path(outdir, "14_run_summary.csv"))
  
  # sessionInfo 저장
  if (isTRUE(cfg$output$save_session_info %||% TRUE)) {
    sink(file.path(outdir, "15_sessionInfo.txt"))
    print(sessionInfo())
    sink()
  }
  
  message("Mfuzz pipeline finished successfully.")
  message("Output directory: ", outdir)
}


# -----------------------------
# [10] auto-run
# -----------------------------
# RStudio에서 Source 버튼만 눌러도 실행되도록 auto-run 처리
args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.yaml")
run_mfuzz_pipeline(config_file)
