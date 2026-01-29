plot_violin_ancova_nice <- function(df, y_var, group_var,
                                    covars = c("subject_age", "PTGENDER", "site", "PTEDUCAT"),
                                    label_style = c("star+p", "star", "pval")) {
  label_style <- match.arg(label_style)
  library(ggplot2)
  library(dplyr)
  library(rlang)
  library(ggpubr)

  # 1) 数据清理并强制 group 为 factor（创建固定 Group 列以稳定 aes）
  df2 <- df %>%
    filter(!is.na(.data[[y_var]]), !is.na(.data[[group_var]])) %>%
    mutate(
      Group = factor(.data[[group_var]], levels = c(1, 2, 3), labels = c("CN", "sMCI", "pMCI"))
    )
  if (nrow(df2) == 0) stop("没有可用数据（y 或 group 全部为 NA）。")
  valid_levels <- levels(df2$Group)
  if (length(valid_levels) < 2) stop("分组少于2个水平，无法比较。")

  # --- Quade 关键步骤（rank + 调整协变量残差） ---
  rank_col <- paste0(y_var, "_rank")
  df2[[rank_col]] <- rank(df2[[y_var]], ties.method = "average")

  covars_present <- covars[covars %in% names(df2)]
  if (length(covars_present) == 0) {
    df2$resid_rank <- df2[[rank_col]]
    quade_regression <- NULL
  } else {
    is_num <- sapply(df2[, covars_present, drop = FALSE], is.numeric)
    num_covs <- covars_present[is_num]
    cat_covs <- covars_present[!is_num]
    rank_covar_names <- character(0)
    if (length(num_covs) > 0) {
      for (cv in num_covs) {
        rn <- paste0(cv, "_rank")
        df2[[rn]] <- rank(df2[[cv]], ties.method = "average")
        rank_covar_names <- c(rank_covar_names, rn)
      }
    }
    rhs_terms <- c(rank_covar_names, cat_covs)
    if (length(rhs_terms) == 0) {
      resid_mod <- lm(as.formula(paste0(rank_col, " ~ 1")), data = df2)
    } else {
      formula_q <- as.formula(paste0(rank_col, " ~ ", paste(rhs_terms, collapse = " + ")))
      resid_mod <- lm(formula_q, data = df2)
    }
    quade_regression <- resid_mod
    df2$resid_rank <- resid(resid_mod)
  }

  # 全局检验与事后
  aov_model <- aov(resid_rank ~ Group, data = df2)
  aov_tab <- anova(aov_model)
  quade_p <- if ("Pr(>F)" %in% colnames(aov_tab)) aov_tab$`Pr(>F)`[1] else NA

  pw <- tryCatch(pairwise.t.test(df2$resid_rank, df2$Group, p.adjust.method = "bonferroni"),
    error = function(e) NULL
  )
  if (is.null(pw)) {
    ph_raw <- data.frame(
      contrast = character(0),
      p.value = numeric(0),
      statistic = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    matp <- pw$p.value
    combs <- which(!is.na(matp), arr.ind = TRUE)
    contrast <- character(0)
    pvals <- numeric(0)
    for (i in seq_len(nrow(combs))) {
      r <- combs[i, 1]
      c <- combs[i, 2]
      g1 <- colnames(matp)[c]
      g2 <- rownames(matp)[r]
      contrast <- c(contrast, paste0(g2, " - ", g1))
      pvals <- c(pvals, matp[r, c])
    }
    ph_raw <- data.frame(contrast = contrast, p.value = pvals, statistic = NA_real_, stringsAsFactors = FALSE)
  }

  # 若没有对比，返回基础图
  if (nrow(ph_raw) == 0) {
    p0 <- ggplot(df2, aes(x = Group, y = .data[[y_var]])) +
      geom_violin(trim = FALSE, fill = "#A6CEE3", width = 0.9) +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "#1F78B4", alpha = 0.8) +
      theme_classic() +
      labs(x = "Group", y = y_var, title = paste0(y_var, " by Group"))
    return(list(plot = p0, quade_regression = quade_regression, aov = aov_model, posthoc = ph_raw, raw_ph = ph_raw))
  }

  # 解析 contrast
  raw_left <- trimws(sapply(strsplit(as.character(ph_raw$contrast), " - "), `[`, 1))
  raw_right <- trimws(sapply(strsplit(as.character(ph_raw$contrast), " - "), `[`, 2))

  fmt_p <- function(p_value) {
    if (is.na(p_value)) return("p=NA")
    if (p_value < 0.001) return("p<0.001")
    paste0("p=", formatC(p_value, digits = 4, format = "f"))
  }

  ph_df <- ph_raw %>%
    mutate(
      group1 = raw_left,
      group2 = raw_right,
      p.value = as.numeric(p.value),
      p.label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      p.text = paste0(p.label, " (", fmt_p(p.value), ")")
    ) %>%
    filter(group1 %in% valid_levels & group2 %in% valid_levels) %>%
    arrange(group1, group2)

  if (nrow(ph_df) == 0) {
    p0 <- ggplot(df2, aes(x = Group, y = .data[[y_var]])) +
      geom_violin(trim = FALSE, fill = "#A6CEE3", width = 0.9) +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "#1F78B4", alpha = 0.8) +
      theme_classic() +
      labs(x = "Group", y = y_var, title = paste0(y_var, " by Group"))
    return(list(plot = p0, quade_regression = quade_regression, aov = aov_model, posthoc = ph_df, raw_ph = ph_raw))
  }

  # xmin/xmax 索引（用于 bracket 横向位置）
  ph_df$xmin <- as.numeric(factor(ph_df$group1, levels = valid_levels))
  ph_df$xmax <- as.numeric(factor(ph_df$group2, levels = valid_levels))

  # 计算每组的高点（95%分位）作为竖线起点（稳健）
  group_top <- tapply(df2[[y_var]], df2$Group, function(x) {
    if (all(is.na(x))) return(NA_real_)
    p <- quantile(x, probs = 0.95, na.rm = TRUE)
    if (is.na(p)) p <- max(x, na.rm = TRUE)
    return(as.numeric(p))
  })
  group_top[is.infinite(group_top) | is.na(group_top)] <- max(df2[[y_var]], na.rm = TRUE)

  # ========== 新堆叠策略：尽量靠近数据上沿，避免拉伸坐标轴 ==========
  y_max_obs <- max(df2[[y_var]], na.rm = TRUE)
  y_min_obs <- min(df2[[y_var]], na.rm = TRUE)
  y_rng <- y_max_obs - y_min_obs
  if (y_rng == 0) y_rng <- abs(y_max_obs) + 1e-6

  # 基线与步长设置为较小的比例，避免上方堆叠过高
  base_above <- y_rng * 0.04
  base_y <- y_max_obs + base_above
  offset <- y_rng * 0.06

  # 按事后顺序线性分配高度
  ph_df <- ph_df %>% mutate(idx = seq_len(n()))
  ph_df$y.position <- base_y + ph_df$idx * offset
  ph_df$label_y <- ph_df$y.position + offset * 0.35
  ph_df$label_x <- (ph_df$xmin + ph_df$xmax) / 2

  if (label_style == "star") {
    ph_df$label_text <- ifelse(ph_df$p.label == "ns", "", ph_df$p.label)
  } else if (label_style == "pval") {
    ph_df$label_text <- fmt_p(ph_df$p.value)
  } else {
    ph_df$label_text <- ph_df$p.text
  }
  ph_df$label_color <- ifelse(ph_df$p.value < 0.05, "red", "black")
  ph_df$label_face <- ifelse(ph_df$p.value < 0.05, "bold", "plain")

  # 绘图：保持小提琴较宽（不要太压缩），并控制坐标轴范围
  color_palette <- c("#E69F00", "#56B4E9", "#009E73")
  p <- ggplot(df2, aes(x = Group, y = .data[[y_var]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.95, width = 0.9) +
    geom_boxplot(
      width = 0.12, outlier.shape = NA,
      fill = c("#E69F00", "#56B4E9", "#009E73"), alpha = 0.6
    ) +
    scale_fill_manual(values = color_palette) +
    scale_x_discrete(limits = valid_levels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(y_min_obs, y_max_obs), clip = "off") +
    theme_classic(base_size = 13) +
    theme(plot.margin = unit(c(1.5, 1, 1, 1), "lines"))

  # Quade 全局 p 用 subtitle
  subtitle_text <- NULL
  subtitle_style <- element_text(size = 10, hjust = 0.5)
  if (!is.na(quade_p)) {
    if (quade_p < 0.001) {
      subtitle_text <- "Quade test p<0.001 ***"
    } else if (quade_p < 0.01) {
      subtitle_text <- paste0("Quade test p=", signif(quade_p, 2), " **")
    } else if (quade_p < 0.05) {
      subtitle_text <- paste0("Quade test p=", signif(quade_p, 2), " *")
    } else {
      subtitle_text <- paste0("Quade test p=", signif(quade_p, 2))
    }
    subtitle_style <- element_text(
      size = 10,
      face = "bold",
      colour = ifelse(quade_p < 0.05, "red", "black"),
      hjust = 0.5
    )
  }

  display_y <- ifelse(tolower(y_var) == "switchrate", "Switch rate", y_var)
  p <- p + labs(
    x = "Group",
    y = display_y,
    fill = "Group",
    title = paste0(display_y, " by groups"),
    subtitle = subtitle_text
  )

  # 画竖线（端点更短但明显）——竖线改细
  vertical_segments1 <- do.call(rbind, lapply(seq_len(nrow(ph_df)), function(i) {
    xmin <- ph_df$xmin[i]
    ytop1 <- group_top[as.character(ph_df$group1[i])] + offset * 0.02
    ybottom <- ph_df$y.position[i] - offset * 0.08
    data.frame(x = xmin, xend = xmin, y = ytop1, yend = ybottom, stringsAsFactors = FALSE)
  }))
  vertical_segments2 <- do.call(rbind, lapply(seq_len(nrow(ph_df)), function(i) {
    xmax <- ph_df$xmax[i]
    ytop2 <- group_top[as.character(ph_df$group2[i])] + offset * 0.02
    ybottom <- ph_df$y.position[i] - offset * 0.08
    data.frame(x = xmax, xend = xmax, y = ytop2, yend = ybottom, stringsAsFactors = FALSE)
  }))

  p <- p + geom_segment(
    data = vertical_segments1, aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, size = 0.6, lineend = "round", color = "black"
  )
  p <- p + geom_segment(
    data = vertical_segments2, aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, size = 0.6, lineend = "round", color = "black"
  )

  # 插入 bracket（横杆）
  p <- p + stat_pvalue_manual(
    ph_df %>% mutate(p.text = ""),
    label = "p.text",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.02,
    bracket.size = 0.9
  )

  # 插入白底 label（更高）
  p <- p + geom_label(
    data = ph_df,
    aes(x = label_x, y = label_y, label = label_text, fontface = label_face, color = label_color),
    inherit.aes = FALSE, fill = "white", size = 4.0,
    label.r = unit(0.16, "lines"), label.padding = unit(0.12, "lines")
  ) +
    scale_color_identity()

  p <- p + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = subtitle_style)

  return(list(plot = p, quade_regression = quade_regression, aov = aov_model, posthoc = ph_df, raw_ph = ph_raw))
}
