# Subsampling-based-feature-selection-for-high-dimensional-data-prediction-
Subsampling Ranking Forward selection (SuRF)—a hybrid framework that combines Elastic-Net penalization, thousands of subsamples for stability ranking, permutation-based p-values, and a greedy forward-selection tree
library(glmnet)
library(tidyverse)
library(caret)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)

# 1.使用 aids 数据集测试数据清洗函数
aids <- read.csv("C:/Users/25478/Downloads/AIDS_clinical_trials.csv", stringsAsFactors = FALSE)
str(aids)
summary(aids)
aids <- aids %>% filter(cd420 > 0)

# 2. 数据清洗与预处理函数
X.c <- NULL
X.o <- aids %>% select(-cd420)
y <- log(aids$cd420 + 1) 

dataclean <- function(X.c, X.o, y) {
  if (!is.null(X.c) && ncol(X.c) > 0) {
    rowSumsXc <- apply(X.c, 1, sum)
    rowSumsXc[rowSumsXc == 0] <- 1
    X.c <- X.c / rowSumsXc
  }
  
  if ((is.null(X.o) || ncol(X.o) == 0) && (!is.null(X.c) && ncol(X.c) > 0)) {
    X <- data.frame(X.c)
  } else if ((is.null(X.o) || ncol(X.o) == 0) && (is.null(X.c) || ncol(X.c) == 0)) {
    stop("No predictor variables given")
  } else if ((is.null(X.c) || ncol(X.c) == 0) && (!is.null(X.o) && ncol(X.o) > 0)) {
    X <- data.frame(X.o)
  } else {
    X <- data.frame(X.c, X.o)
  }
  
  Data.Xy <- data.frame(X, y)
  return(Data.Xy)
}

# 2.1 Subsample.w：一次子抽样及 LASSO 模型拟合
Subsample.w <- function(data, fold, Alpha, prop, weights, family, Type) {
  n <- nrow(data)
  p <- ncol(data)
  ID <- seq_len(n)
  
  if (family$family == "binomial") {
    if (is.null(weights) || identical(weights, FALSE)) {
      w <- rep(1, n)
    } else {
      N1 <- sum(data[, p] == 1)
      N0 <- sum(data[, p] == 0)
      w <- ifelse(data[, p] == 1, 1, N1 / N0)
    }
  } else {
    w <- if (is.null(weights) || identical(weights, FALSE)) rep(1, n) else weights
  }
  
  test.ID <- sample(ID, size = ceiling(n * prop))
  train.ID <- setdiff(ID, test.ID)
  
  train.X <- data[train.ID, -p, drop = FALSE]
  train.Y <- data[train.ID, p]
  test.X  <- data[test.ID, -p, drop = FALSE]
  test.Y  <- data[test.ID, p]
  
  if (family$family == "binomial") {
    train.Y <- as.factor(train.Y)
    test.Y  <- as.factor(test.Y)
  }
  
  NonZero_ind <- which(colSums(train.X) != 0)
  
  if (is.null(weights) || identical(weights, FALSE)) {
    model.cv <- cv.glmnet(as.matrix(train.X[, NonZero_ind, drop = FALSE]), train.Y,
                          family = family$family, type.measure = Type, nfolds = fold,
                          alpha = Alpha, standardize = TRUE, parallel = FALSE)
  } else {
    model.cv <- cv.glmnet(as.matrix(train.X[, NonZero_ind, drop = FALSE]), train.Y,
                          family = family$family, type.measure = Type, nfolds = fold,
                          alpha = Alpha, standardize = TRUE, weights = w[train.ID], parallel = FALSE)
  }
  
  lambda_min <- model.cv$lambda.min
  preds <- predict(model.cv, as.matrix(test.X[, NonZero_ind, drop = FALSE]), s = "lambda.min")
  
  if (family$family == "binomial") {
    pred_class <- as.numeric(preds > 0.5)
    true_class <- as.numeric(test.Y)
    error_rate <- mean(pred_class != true_class)
  } else {
    error_rate <- mean((as.numeric(preds) - test.Y)^2)
  }
  
  coef_sparse <- coef.glmnet(model.cv, s = lambda_min)
  coef_mat <- as.matrix(coef_sparse)
  coef_vals <- coef_mat[, 1]
  nonzero_idx <- which(coef_vals != 0)
  COEF <- data.frame(Variable = rownames(coef_mat)[nonzero_idx],
                     Coefficient = coef_vals[nonzero_idx],
                     stringsAsFactors = FALSE)
  
  beta_vec <- rep(0, ncol(train.X))
  beta_vec[NonZero_ind] <- as.vector(coef_sparse[-1])
  beta_vec <- c(beta_vec, as.numeric(coef_sparse[1]))
  
  list(lambda = lambda_min, coef = COEF, Error = error_rate, Beta = beta_vec)
}

# 2.2 Subsample_B：重复 B 次子抽样（支持并行）
Subsample_B <- function(B, data, fold, Alpha, prop, weights, ncores, family) {
  Type <- if (family$family == "binomial") "class" else "deviance"
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    clusterExport(cl, varlist = c("Subsample.w"), envir = environment())
  }
  p <- ncol(data) - 1
  RESULTS <- foreach::foreach(i = 1:B, .combine = rbind, .packages = "glmnet") %dopar% {
    model <- Subsample.w(data, fold, Alpha, prop, weights, family, Type)
    c(model$Error, model$lambda, model$Beta)
  }
  if (ncores > 1) stopCluster(cl)
  avg_error <- mean(RESULTS[, 1])
  lambda_vec <- RESULTS[, 2]
  beta_matrix <- RESULTS[, -c(1, 2), drop = FALSE]
  list(Class.Err = avg_error, Lambda = lambda_vec, BETA = beta_matrix)
}


# 2.3 Ranking：对变量进行排序
Ranking <- function(data, model) {
  BETA <- model$BETA
  binary_beta <- (BETA != 0) * 1
  binary_predictors <- binary_beta[, -ncol(binary_beta), drop = FALSE]
  freq <- colSums(binary_predictors)
  predictor_names <- colnames(data)[-ncol(data)]
  ranking_order <- order(freq, decreasing = TRUE)
  rank_table <- data.frame(variable = predictor_names[ranking_order],
                           frequency = freq[ranking_order],
                           stringsAsFactors = FALSE)
  rank_table <- subset(rank_table, frequency > 0)
  
  if (nrow(rank_table) > 1) {
    X_sel <- data[, -ncol(data), drop = FALSE][, rank_table$variable, drop = FALSE]
    corr_mat <- stats::cor(X_sel)
    to_remove <- NULL
    for (i in seq_len(ncol(corr_mat))) {
      if (i < ncol(corr_mat)) {  # 避免下标越界
        high_corr <- which(corr_mat[i, (i+1):ncol(corr_mat)] > 0.999)
        if (length(high_corr) > 0) {
          to_remove <- c(to_remove, i + high_corr)
        }
      }
    }
    if (length(to_remove) > 0) {
      keep_idx <- setdiff(seq_len(nrow(rank_table)), unique(to_remove))
      rank_table <- rank_table[keep_idx, , drop = FALSE]
    }
  }
  return(rank_table)
}



# 2.4 update_dev：置换检验计算 deviance 分布
update_dev <- function(data, vslist, C, weights, ncores, family) {
  p <- ncol(data)
  N <- nrow(data)
  
  if (family$family == "binomial") {
    if (is.null(weights) || identical(weights, FALSE)) {
      w <- rep(1, N)
    } else {
      N1 <- sum(data[, p] == 1)
      N0 <- sum(data[, p] == 0)
      w <- ifelse(data[, p] == 1, 1, N1 / N0)
    }
  } else {
    w <- if (is.null(weights) || identical(weights, FALSE)) rep(1, N) else weights
  }
  
  X <- data[, -p, drop = FALSE]
  y <- data[, p]
  
  if (length(vslist) > 0) {
    xmat2 <- data.frame(X[, vslist, drop = FALSE], y = y)
  } else {
    xmat2 <- data.frame(y = y)
  }
  mod2 <- glm(y ~ ., data = xmat2, family = family, weights = w, control = list(maxit = 500))
  
  dev <- rep(0, C)
  nonzeroind <- which(colSums(X) != 0)
  
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    dev <- foreach::foreach(j = 1:C, .combine = c, .packages = "stats") %dopar% {
      ind <- sample(seq_len(N))
      XX <- X
      yy <- y
      if (length(vslist) > 0) {
        XX[, vslist] <- X[ind, vslist]
      }
      yy <- y[ind]
      ww <- w[ind]
      devseq <- sapply(setdiff(nonzeroind, vslist), function(col_idx) {
        xvar <- XX[, col_idx]
        dat1 <- data.frame(XX[, vslist, drop = FALSE], xvar = xvar, y = yy)
        mod1 <- glm(y ~ ., data = dat1, family = family, weights = ww, control = list(maxit = 500))
        return(deviance(mod2) - deviance(mod1))
      })
      max(devseq, na.rm = TRUE)
    }
    stopCluster(cl)
  } else {
    for(j in 1:C){
      ind <- sample(seq_len(N))
      XX <- X
      yy <- y
      if (length(vslist) > 0) {
        XX[, vslist] <- X[ind, vslist]
      }
      yy <- y[ind]
      ww <- w[ind]
      devseq <- sapply(setdiff(nonzeroind, vslist), function(col_idx) {
        xvar <- XX[, col_idx]
        dat1 <- data.frame(XX[, vslist, drop = FALSE], xvar = xvar, y = yy)
        mod1 <- glm(y ~ ., data = dat1, family = family, weights = ww, control = list(maxit = 500))
        return(deviance(mod2) - deviance(mod1))
      })
      dev[j] <- max(devseq, na.rm = TRUE)
    }
  }
  return(dev)
}

# 3.1 selectnew：从候选变量中筛选新变量
selectnew <- function(vslist, ranktable, data, weights, ncores = 1, family, alpha_l, alpha_u, C) {
  ecdf_fun <- function(x, perc) stats::ecdf(x)(perc)
  n <- nrow(data)
  p <- ncol(data)
  
  if (family$family == "binomial") {
    if (is.null(weights) || identical(weights, FALSE)) {
      w <- rep(1, n)
    } else {
      N1 <- sum(data[, p] == 1)
      N0 <- sum(data[, p] == 0)
      w <- ifelse(data[, p] == 1, 1, N1 / N0)
    }
  } else {
    w <- if (is.null(weights) || identical(weights, FALSE)) rep(1, n) else weights
  }
  
  new_vars <- setdiff(ranktable$variable, vslist)
  c.ranktable <- ranktable[ranktable$variable %in% new_vars, ]
  vclist <- as.character(c.ranktable$variable)
  freq <- c.ranktable$frequency
  
  # 避免嵌套并行：将 update_dev 的 ncores 固定为 1
  newdev.dist <- update_dev(data, vslist, C, weights, ncores = 1, family)
  
  if (is.null(vslist)) {
    vtlist <- matrix(vclist, ncol = 1)
  } else {
    vtlist <- t(sapply(vclist, function(x) c(vslist, x)))
  }
  colnames(vtlist) <- paste("V", 1:ncol(vtlist), sep = "")
  
  pval_results <- apply(vtlist, 1, function(candidate_vars) {
    dat <- data[, c(candidate_vars, names(data)[p]), drop = FALSE]
    colnames(dat)[ncol(dat)] <- "y"
    mod <- glm(y ~ ., data = dat, family = family, weights = w, control = list(maxit = 500))
    vardev <- deviance(mod)
    r <- sum(newdev.dist >= vardev)
    # 2) 用 (r+1)/(C+1) 平滑
    pval <- (r + 1) / (length(newdev.dist) + 1)
    c(pval, vardev)
  })
  pval_table <- data.frame(vtlist, freq = freq, t(pval_results))
  colnames(pval_table)[(ncol(pval_table)-1):ncol(pval_table)] <- c("pval", "deviance")
  
  pval_table <- pval_table[order(-pval_table$freq, -pval_table$deviance), ]
  cutoff <- quantile(newdev.dist, 1 - alpha_u)
  pval_table <- pval_table[pval_table$deviance > cutoff, ]
  
  return(pval_table)
}


# 3.2 selpath：递归构建变量选择路径
selpath <- function(data, weights, ranktable, ncores, family, C, alpha_u, max_depth = 4) {
  # 初始化：第一层节点
  current_depth <- 1
  all_nodes <- list()       # 存储所有生成的节点
  current_nodes <- list(list(alpha.range = c(0, alpha_u), selvar = NULL, vslist = NULL))
  all_nodes <- current_nodes
  
  while (length(current_nodes) > 0 && current_depth < max_depth) {
    next_nodes <- list()
    cat("Current recursion depth:", current_depth, "The current number of nodes:", length(current_nodes), "\n")
    for (node in current_nodes) {
      vslist <- node$vslist
      # 检查剩余候选变量
      remaining <- setdiff(ranktable$variable, vslist)
      if (length(remaining) == 0) {
        next
      } else {
        # 调用 selectnew 获取当前节点扩展的候选变量（注意：这里 ncores 固定为 1，避免嵌套并行）
        out <- selectnew(vslist = vslist, ranktable = ranktable, data = data, 
                         weights = weights, ncores = 1, family = family, 
                         alpha_l = node$alpha.range[1], alpha_u = node$alpha.range[2], C = C)
        if (nrow(out) > 0) {
          for (j in seq_len(nrow(out))) {
            candidate_var <- if (is.null(vslist)) out[j, "V1"] else out[j, paste0("V", length(vslist) + 1)]
            new_vslist <- if (is.null(vslist)) candidate_var else c(vslist, candidate_var)
            new_node <- list(alpha.range = node$alpha.range,
                             selvar = candidate_var,
                             vslist = new_vslist,
                             pval = out[j, "pval"],
                             dev = out[j, "deviance"])
            next_nodes[[length(next_nodes) + 1]] <- new_node
          }
        }
      }
    }
    cat("The number of nodes added to this layer:", length(next_nodes), "\n")
    all_nodes <- c(all_nodes, next_nodes)
    current_nodes <- next_nodes
    current_depth <- current_depth + 1
  }
  # 返回最后一层生成的节点作为选择路径的结果
  return(list(selpoint = current_nodes, sel.nodes = all_nodes))
}


# 3.3 selvar_alpha：选择最终变量集合
selvar_alpha <- function(res, alpha) {
  valid_nodes <- lapply(res$selpoint, function(node) {
    if (!is.null(node$alpha.range) && node$alpha.range[1] <= alpha && node$alpha.range[2] >= alpha)
      return(node)
    else
      return(NULL)
  })
  valid_nodes <- Filter(Negate(is.null), valid_nodes)
  if (length(valid_nodes) == 0) {
    return(list(vslist = NULL, alpha.range = c(alpha, alpha)))
  }
  final_node <- valid_nodes[[which.max(sapply(valid_nodes, function(x) length(x$vslist)))]]
  return(final_node)
}


# 4. SURF 主函数
SURF <- function(Xo, y, X = NULL, fold = 10, Alpha = 0.5, prop = 0.1, weights = FALSE,
                 B = 1000, C = 400, ncores = ncores_used, display.progress = TRUE,
                 family = stats::gaussian(link = "identity"), alpha_u = 0.5, alpha = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
)
  {
  traindata <- dataclean(X.c = X, X.o = Xo, y = y)
  if (display.progress) print("Data cleaning completed.")
  
  tempmod <- Subsample_B(B = B, data = traindata, fold = fold, Alpha = Alpha, prop = prop,
                         weights = weights, ncores = ncores, family = family)
  if (display.progress) print("Subsampling completed.")
  
  rank_table <- Ranking(traindata, tempmod)
  if (display.progress) print("Ranking completed.")
  
  if (C == 0) {
    ans <- list(Bmod = tempmod, trdata = traindata, rank_table = rank_table, family = family)
    class(ans) <- "SURF_fit"
    warning("SURF only performs ranking without permutation-based selection.")
    return(ans)
  } else {
    modpath <- selpath(data = traindata, weights = weights, ranktable = rank_table,
                       ncores = ncores, family = family, C = C, alpha_u = alpha_u)
    if (display.progress) print("Selection path completed.")
    
    valid_alpha <- alpha[alpha > 0 & alpha <= alpha_u]
    if (length(valid_alpha) == 0) {
      print("Warning: no valid alpha provided; using alpha_u")
      valid_alpha <- alpha_u
    }
    selmod <- lapply(valid_alpha, function(a) selvar_alpha(res = modpath, alpha = a))
    for (i in seq_along(valid_alpha)) {
      selmod[[i]]$alpha <- valid_alpha[i]
    }
    if (display.progress) print("Final variable selection completed.")
    
    ans <- list(Bmod = tempmod, trdata = traindata, rank_table = rank_table,
                modpath = modpath, selmod = selmod, family = family)
    class(ans) <- "SURF_fit"
    return(ans)
  }
}

ncores_used <- 2
cl <- makeCluster(ncores_used)
clusterEvalQ(cl, {
  library(glmnet)
  library(foreach)
})
clusterExport(cl, varlist = c("Subsample.w", "update_dev", "selectnew", "Ranking",
                              "selpath", "selvar_alpha", "dataclean"), envir = environment())
registerDoParallel(cl)

# 调用 SURF 
surf_result <- SURF(Xo = X.o, y = y, X = NULL,
                    B = 1000,
                    C = 400,
                    ncores = ncores_used,
                    display.progress = TRUE,
                    family = stats::gaussian(link = "identity"),
                    alpha_u = 0.5,
                    alpha = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
)

stopCluster(cl)
cat("Variable set selected by SURF (at each alpha level) ：\n")
print(surf_result$selmod)

# 5. 提取最终变量并构建最终预测模型
final_vars <- surf_result$selmod[[1]]$vslist
cat("Final selected var ：\n")
print(final_vars)
train_data <- surf_result$trdata
colnames(train_data)[ncol(train_data)] <- "y"

# 6. 构造公式并拟合线性模型
formula_str <- paste("y ~", paste(final_vars, collapse = " + "))
final_formula <- as.formula(formula_str)
final_lm <- lm(final_formula, data = train_data)
summary(final_lm)

# 7. 在同一数据上做预测并计算 MSE
preds <- predict(final_lm, newdata = train_data)
mse   <- mean((preds - train_data$y)^2)
cat("Training set MSE:", mse, "\n")
rmse <- sqrt(mse)
cat("Training set RMSE:", rmse, "\n")
r2 <- summary(final_lm)$r.squared
cat("Training set R-squared:", r2, "\n")

plot_data <- data.frame(
  Fitted    = predict(final_lm, newdata = train_data),
  Residuals = train_data$y - predict(final_lm, newdata = train_data)
)

# Residuals vs Fitted
ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs Fitted",
    x     = "Fitted log(cd420 + 1)",
    y     = "Residuals"
  ) +
  theme_minimal()

# QQ-plot 检查残差正态性
qqnorm(residuals(final_lm))
qqline(residuals(final_lm), col = "red", lwd = 2)




