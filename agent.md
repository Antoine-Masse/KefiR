# Ce package KefiR doit être publié sur gitgub (nouvelle mise à jour à venir).

Il s'agit pour toi de modifier le moins possible les fichiers, d'en créer le moins possible.
De conserver au maximum le style, la doc....

Ta seul tache est de faire tous les contrôles roxygen2 nécessaires et modifications minimales pour que puisse avoir lieu le git (fait par l'utilisateur).


Exemple du dernier bilan roxygen2 fait manuellement : 
  <environment: 0x000001fbfeb40900>
  [1] "Making evolutiv approach on  4000  iterations."
  De :  201 - 400 à : 401 - 600 -  0.1381264 pour  6134  parents en prélèvements  53 .
  De :  201 - 400 à : 601 - 800 -  1.922794e-08 pour  2621  parents en prélèvements  53 .
  De :  601 - 800 à : 801 - 1000 -  3.39926e-07 pour  1856  parents en prélèvements  10 .
  De :  801 - 1000 à : 1001 - 1200 -  0.8135395 pour  1297  parents en prélèvements  8 .
  De :  801 - 1000 à : 1201 - 1400 -  0.01015226 pour  778  parents en prélèvements  8 .
  [1] "Evolutive approach really made on  2138 iterations."
  [1] "il reste parents :  4"
  [1] "Meilleur modèle sélectionné :"
  mpg ~ I(1/hp) + I(1/wt)
  <environment: 0x000001fbfeb40900>
  Capacité de prédiction brute :  0.88 
  Capacité de prédiction médiane (bootstrap):  0.89 
  Capacité de prédiction inférieure 95% :  0.85 
  > myreg2 <- evolreg(mtcars,"cyl")
  Error in ks.test.default(x_i, x_j) : not enough 'y' data
  Calls: evolreg -> dvar -> m.test -> .one_factor_analysis -> pairwise
  Execution halted

❯ checking code files for non-ASCII characters ... WARNING
  Found the following files with non-ASCII characters:
    R/catego.R
    R/cooks.distance_lmer.R
    R/evolreg.r
    R/gage_rr.R
    R/m.test.R
    R/pairwise.R
    R/pairwise.boot.R
    R/pde.R
    R/sample_bootstrap.R
    R/sys_analyze_ancova_structure.R
    R/sys_ancova_analysis.R
    R/sys_auto_preprocess_g.R
    R/sys_control_independence.R
    R/sys_detect_model_type.R
    R/sys_filter_small_groups.R
    R/sys_formulator.R
    R/sys_formulator_safe.R
    R/sys_manova_analysis.R
    R/sys_mixed_model_analysis.R
    R/sys_multi_factor_analysis.R
    R/sys_normality.R
    R/sys_one_factor_analysis.R
    R/sys_pairwise_medpb2.R
    R/sys_plot_with_letters.R
    R/sys_posthoc.R
    R/sys_posthoc_ANCOVA.R
    R/sys_posthoc_MANOVA.R
    R/sys_posthoc_mixed_model.R
    R/sys_select_ss_type.R
    R/sys_variance.R
    R/valreg.R
  Portable packages must use only ASCII characters in their R code and
  NAMESPACE directives, except perhaps in comments.
  Use \uxxxx escapes for other characters.
  Function 'tools::showNonASCIIfile' can help in finding non-ASCII
  characters in files.

❯ checking for code/documentation mismatches ... WARNING
  Codoc mismatches from Rd file 'm.test.Rd':
  m.test
    Code : function(x = NULL, g = NULL, data = NULL, formula = NULL,
                   paired = FALSE, id = NULL, wt = NULL, within = NULL,
                   between = NULL, alpha = 0.05, control = NULL, verbose
                   = TRUE, plot = TRUE, return = TRUE, boot = TRUE,
                   boot_type = NULL, iter = 0, conf = 0.95, maxcat = 50,
                   silent = TRUE, code = FALSE, debug = FALSE)
    Docs : function(formula, data = NULL, alpha = 0.05, verbose = TRUE,
                   return = TRUE, id = NULL, control = c(), maxcat = 50,
                   plot = TRUE, silent = TRUE, boot = TRUE, iter = 500,
                   conf = 0.95, code = FALSE, debug = FALSE)
    Noms d'arguments dans le code mais pas dans la doc :
      x g paired wt within between boot_type
    Incohérence dans les noms d'arguments (les 3 premiers) :
      Position: 1 Code: x Docs: formula
      Position: 2 Code: g Docs: data
      Position: 3 Code: data Docs: alpha
    Incohérence dans les valeurs par défaut des arguments :
      Name: 'formula' Code: NULL Docs: 
      Name: 'control' Code: NULL Docs: c()
      Name: 'iter' Code: 0 Docs: 500
  
  Codoc mismatches from Rd file 'pairwise.Rd':
  pairwise
    Code : function(x, g, type = "mean", alpha = 0.05, control = c(),
                   pool.sd = FALSE, silent = TRUE, boot = FALSE,
                   boot_type = mean, tr = 0.2, iter = 500, conf = 0.95,
                   paired = FALSE, debug = FALSE)
    Docs : function(x, g, type = "mean", alpha = 0.05, control = c(),
                   pool.sd = FALSE, silent = TRUE, boot = FALSE, iter =
                   500, conf = 0.95, debug = FALSE, p.adjust.method =
                   "holm", conf.level = 0.95, alternative = "two.sided")
    Noms d'arguments dans le code mais pas dans la doc :
      boot_type tr paired
    Noms d'arguments dans la doc mais pas dans le code :
      p.adjust.method conf.level alternative
    Incohérence dans les noms d'arguments (les 3 premiers) :
      Position: 9 Code: boot_type Docs: iter
      Position: 10 Code: tr Docs: conf
      Position: 11 Code: iter Docs: debug
  
  Codoc mismatches from Rd file 'pairwise.boot.Rd':
  pairwise.boot
    Code : function(x, g, iter = 500, mu = "mean", paired = FALSE, debug
                   = FALSE, verbose = F, return = T)
    Docs : function(x, g, iter = 500, mu = "meanbp", debug = FALSE,
                   p.adjust.method = "holm", conf.level = 0.95,
                   alternative = "two.sided")
    Noms d'arguments dans le code mais pas dans la doc :
      paired verbose return
    Noms d'arguments dans la doc mais pas dans le code :
      p.adjust.method conf.level alternative
    Incohérence dans les noms d'arguments (les 3 premiers) :
      Position: 5 Code: paired Docs: debug
      Position: 6 Code: debug Docs: p.adjust.method
      Position: 7 Code: verbose Docs: conf.level
    Incohérence dans les valeurs par défaut des arguments :
      Name: 'mu' Code: "mean" Docs: "meanbp"
  
  Codoc mismatches from Rd file 'valreg.Rd':
  valreg
    Code : function(reg, verbose = TRUE, nvar = NULL, boot = TRUE, alpha
                   = 0.05, conf.level = 0.95, plot = FALSE, data = c(),
                   raintest_alpha = 0.05, dwtest_alpha = 0.03,
                   shapiro_alpha = 0.05, bptest_alpha = 0.05, k = 0,
                   orderDW = NULL, tolerance = "basic", debug = FALSE)
    Docs : function(reg, verbose = TRUE, nvar = 5, boot = TRUE, alpha =
                   0.05, conf.level = 0.95, plot = FALSE, data = c(),
                   code = FALSE, raintest_alpha = 0.05, dwtest_alpha =
                   0.03, shapiro_alpha = 0.05, bptest_alpha = 0.05)
    Noms d'arguments dans le code mais pas dans la doc :
      k orderDW tolerance debug
    Noms d'arguments dans la doc mais pas dans le code :
      code
    Incohérence dans les noms d'arguments (les 3 premiers) :
      Position: 9 Code: raintest_alpha Docs: code
      Position: 10 Code: dwtest_alpha Docs: raintest_alpha
      Position: 11 Code: shapiro_alpha Docs: dwtest_alpha
    Incohérence dans les valeurs par défaut des arguments :
      Name: 'nvar' Code: NULL Docs: 5

❯ checking files in 'vignettes' ... WARNING
  Files in the 'vignettes' directory but no files in 'inst/doc':
    'comment-utiliser-mon-package.Rmd'
  Files named as vignettes but with no recognized vignette engine:
     'vignettes/comment-utiliser-mon-package.Rmd'
  (Is a VignetteBuilder field missing?)

❯ checking package dependencies ... NOTE
  Imports includes 30 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking for hidden files and directories ... NOTE
  Found the following hidden files and directories:
    .testeur.R
    inst/dev/R_poubelle/.align_pairs.R
    inst/dev/R_poubelle/.auto_preprocess_g.R
    inst/dev/R_poubelle/.boots.R
    inst/dev/R_poubelle/.control_independence.R
    inst/dev/R_poubelle/.dbg.R
    inst/dev/R_poubelle/.detect_model_type.R
    inst/dev/R_poubelle/.drop_error_term.R
    inst/dev/R_poubelle/.exit.R
    inst/dev/R_poubelle/.format_pval.R
    inst/dev/R_poubelle/.formulator.R
    inst/dev/R_poubelle/.formulator_safe.R
    inst/dev/R_poubelle/.manova_analysis.R
    inst/dev/R_poubelle/.msg.R
    inst/dev/R_poubelle/.multi_factor_analysis02.R
    inst/dev/R_poubelle/.normality.R
    inst/dev/R_poubelle/.normalize_formula_dollar.R
    inst/dev/R_poubelle/.normalize_from_formula.R
    inst/dev/R_poubelle/.one_factor_analysis.R
    inst/dev/R_poubelle/.posthoc.R
    inst/dev/R_poubelle/.strip_data_dollar_safe.R
    inst/dev/R_poubelle/.testeur.R
    inst/dev/R_poubelle/.variance.R
    inst/dev/R_poubelle/.vbse.R
    inst/dev/Save/.dbg.zip
    inst/dev/Save/.detect_model_type.zip
    vignettes/.RData
  These were most likely included in error. See section 'Package
  structure' in the 'Writing R Extensions' manual.

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking DESCRIPTION meta-information ... NOTE
  Packages listed in more than one of Depends, Imports, Suggests, Enhances:
    'WRS2' 'tseries'
  A package should be listed in only one of these fields.

❯ checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    'INSTALLATION_GUIDE.md'

❯ checking dependencies in R code ... NOTE
  importations '::' ou ':::' non déclarées depuis :
    'broom' 'dae' 'rcompanion' 'rlang' 'rstatix' 'withr'
  appels 'loadNamespace' ou 'requireNamespace' non déclarés depuis :
    'dae' 'multcomp' 'rcompanion' 'rlang' 'rstatix' 'withr'
  Objets manquants ou non exportés :
    'emmeans::cld' 'lmerTest::anova'
  appel ':::' qui devrait être '::' : 'lmerTest:::as_lmerModLmerTest'
    See the note in ?`:::` about the use of this operator.
  Objet non exporté importé par appel ':::' : 'lmerTest:::summary.lmerModLmerTest'
    See the note in ?`:::` about the use of this operator.

❯ checking R code for possible problems ... [74s] NOTE
  .analyze_ancova_structure: no visible global function definition for
    'terms'
  .ancova_analysis: no visible global function definition for 'residuals'
  .boots: no visible global function definition for 'median'
  .boots: no visible global function definition for 'quantile'
  .detect_model_type: no visible global function definition for 'terms'
  .diagnostic_influence: no visible global function definition for 'coef'
  .diagnostic_influence: no visible global function definition for
    'hatvalues'
  .diagnostic_influence: no visible global function definition for
    'dfbetas'
  .drop_error_term: no visible global function definition for 'terms'
  .drop_error_term: no visible global function definition for
    'delete.response'
  .drop_error_term: no visible global function definition for
    'reformulate'
  .manova_analysis: no visible global function definition for 'setNames'
  .manova_analysis: no visible global function definition for
    'complete.cases'
  .manova_analysis: no visible global function definition for 'terms'
  .manova_analysis: no visible global function definition for 'cov'
  .manova_analysis: no visible global function definition for
    'mahalanobis'
  .manova_analysis: no visible global function definition for 'qchisq'
  .manova_analysis: no visible global function definition for 'manova'
  .mixed_model_analysis: no visible global function definition for
    'terms'
  .mixed_model_analysis: possible error in valreg(reg = mixed_model,
    verbose = verbose, code = code, alpha = alpha, boot = FALSE, plot =
    FALSE, k = k, tolerance = "extrem", orderDW = NULL, debug = debug):
    argument inutilisé (code = code)
  .multi_factor_analysis : get_residuals: no visible global function
    definition for 'residuals'
  .multi_factor_analysis : get_residuals: no visible binding for global
    variable 'residuals'
  .multi_factor_analysis : get_studentized_residuals: no visible global
    function definition for 'rstudent'
  .multi_factor_analysis : get_studentized_residuals: no visible global
    function definition for 'rstandard'
  .multi_factor_analysis : get_studentized_residuals : <anonymous>: no
    visible global function definition for 'residuals'
  .multi_factor_analysis: no visible global function definition for
    'terms'
  .multi_factor_analysis: no visible global function definition for 'BIC'
  .multi_factor_analysis: no visible global function definition for
    'reshape'
  .multi_factor_analysis: no visible global function definition for
    'friedman.test'
  .multi_factor_analysis: no visible global function definition for
    'residuals'
  .normality: no visible global function definition for 'p.adjust'
  .normalize_from_formula: no visible binding for global variable
    'na.pass'
  .normalize_from_formula: no visible global function definition for
    'complete.cases'
  .one_factor_analysis : .safe_identify_outliers: no visible global
    function definition for 'quantile'
  .one_factor_analysis: no visible global function definition for
    'mood.test'
  .one_factor_analysis: no visible global function definition for
    'ansari.test'
  .one_factor_analysis: no visible binding for global variable 'var'
  .pairwise_medpb2: no visible global function definition for 'p.adjust'
  .plot_with_letters: no visible global function definition for 'boxplot'
  .plot_with_letters: no visible global function definition for
    'stripchart'
  .posthoc: possible error in .dbg(".posthoc() - Student with control.",
    ".posthoc() - Student avec présence d'un Témoin.", debug = debug,
    verbose = verbose): argument inutilisé (verbose = verbose)
  .posthoc: no visible global function definition for 'str_split_fixed'
  .posthoc: no visible global function definition for 'games_howell_test'
  .posthoc : mgg_test: no visible global function definition for 'median'
  .posthoc: no visible global function definition for 'binom.test'
  .posthoc: no visible global function definition for 'median'
  .posthoc : .safe_identify_outliers: no visible global function
    definition for 'quantile'
  .posthoc : outlier: no visible global function definition for 'median'
  .posthoc: no visible global function definition for 'quantile'
  .posthoc_ANCOVA: no visible global function definition for 'terms'
  .posthoc_ANCOVA: no visible global function definition for 'formula'
  .select_ss_type: no visible global function definition for 'terms'
  .variance: no visible global function definition for 'formula'
  Mode: no visible global function definition for 'density'
  auto_preprocess_g: no visible global function definition for 'exit'
  auto_preprocess_g: no visible global function definition for 'msg'
  auto_preprocess_g : create_bins: no visible global function definition
    for 'quantile'
  auto_preprocess_g: no visible binding for global variable 'verbose'
  auto_preprocess_g: no visible binding for global variable 'k'
  biplt: no visible global function definition for 'rainbow'
  biplt: no visible global function definition for 'grid'
  bootreg: no visible global function definition for 'getCall'
  bootreg: no visible global function definition for 'formula'
  bootreg: no visible global function definition for 'get_all_vars'
  bootreg : mode: no visible global function definition for 'density'
  bootreg: no visible binding for global variable 'median'
  bootreg : boxplot_Pr: no visible binding for global variable 'quantile'
  bootreg : boxplot_Pr: no visible global function definition for
    'boxplot'
  bootreg : boxplot_Pr: no visible global function definition for
    'abline'
  bootreg : <anonymous>: no visible global function definition for
    'median'
  bootreg: no visible global function definition for 'boxplot'
  check_ech: no visible global function definition for 'rnorm'
  control_independence: no visible global function definition for 'exit'
  control_independence: no visible global function definition for
    'p.adjust'
  cooks.distance_lmer: no visible global function definition for
    'influence'
  cor.e: no visible global function definition for 'rnorm'
  cor_ech: no visible global function definition for 'quantile'
  corrigraph: no visible binding for global variable 'var'
  corrigraph : check_mdl: no visible binding for global variable 'var'
  corrigraph : check_mdl: no visible global function definition for 'BIC'
  corrigraph: no visible global function definition for 'BIC'
  corrigraph: no visible binding for global variable 'weight'
  corrigraph: no visible global function definition for 'quantile'
  corrigraph: no visible global function definition for 'var'
  corrigraph: no visible global function definition for 'as.dist'
  corrigraph: no visible global function definition for 'hclust'
  corrigraph: no visible global function definition for 'cutree'
  dsc: no visible global function definition for 'get_all_vars'
  dsc: no visible global function definition for 'formula'
  dsc : <anonymous>: no visible global function definition for 'rbeta'
  dsc2: no visible global function definition for 'formula'
  dsc2: no visible global function definition for 'get_all_vars'
  dsc2: no visible global function definition for 'error'
  dsc2: no visible global function definition for 'boxplot'
  dvar: no visible global function definition for 'chisq.test'
  dvar: no visible global function definition for 'qbinom'
  dvar: no visible global function definition for 'get_all_vars'
  evolreg: no visible global function definition for 'formula'
  evolreg: no visible global function definition for 'get_all_vars'
  evolreg: no visible global function definition for 'BIC'
  evolreg: no visible global function definition for 'glm'
  evolreg: no visible global function definition for 'binomial'
  evolreg: no visible binding for global variable 'logit'
  evolreg: no visible global function definition for 'drop1'
  evolreg: no visible global function definition for 'quantile'
  evolreg: no visible global function definition for 'median'
  evolreg: no visible binding for global variable 'form'
  exp_dyn: no visible global function definition for 'median'
  identify_ech: no visible global function definition for 'rnorm'
  int.ech: no visible global function definition for 'var'
  int.ech: no visible global function definition for 'qt'
  int.pop: no visible global function definition for 'qnorm'
  int.prop: no visible global function definition for 'qnorm'
  int.prop.table: warning in matrix(myIC, nc = ncol(x), nr = nrow(x)):
    correspondance partielle d'argument de 'nr' par rapport à 'nrow'
  int.prop.table: warning in matrix(myIC, nc = ncol(x), nr = nrow(x)):
    correspondance partielle d'argument de 'nc' par rapport à 'ncol'
  int.prop.table: warning in matrix(myIC, nc = ncol(x), nr = nrow(x),
    byrow = TRUE): correspondance partielle d'argument de 'nr' par
    rapport à 'nrow'
  int.prop.table: warning in matrix(myIC, nc = ncol(x), nr = nrow(x),
    byrow = TRUE): correspondance partielle d'argument de 'nc' par
    rapport à 'ncol'
  jb.norm.test: no visible global function definition for 'rnorm'
  jb.normtest: no visible global function definition for 'rnorm'
  kurtosis.norm.test: no visible global function definition for 'rnorm'
  lms.to.table: no visible binding for global variable 'reg'
  m.test: no visible global function definition for 'terms'
  m.test: no visible binding for global variable 'na.pass'
  m.test: no visible global function definition for 'model.response'
  m.test: no visible global function definition for 'update'
  m.test: no visible global function definition for 'complete.cases'
  m.test: possible error in .posthoc_mixed_model(mixed_model =
    mixed_model, significant_effects = sig_effects, alpha = alpha, method
    = "tukey", conf.level = conf, verbose = verbose, code = code, debug =
    debug, code = code, k = k): argument formel "code" correspondant à
    plusieurs arguments fournis
  m.test: no visible binding for global variable 'pvals'
  m.test: no visible binding for global variable 'median'
  m.test: no visible binding for global variable 'IQR'
  meanbp: no visible global function definition for 'quantile'
  meanbp: no visible global function definition for 'median'
  mm.test: no visible binding for global variable 'i'
  mm.test: no visible global function definition for 'pairwise.t.test'
  mm.test: no visible global function definition for
    'pairwise.wilcox.test'
  pairwise : lincon_to_pairwise: warning in matrix(rep(NA,
    (length(t1$fnames) - 1)^2), nc = length(t1$fnames) - 1, nr =
    length(t1$fnames) - 1): correspondance partielle d'argument de 'nr'
    par rapport à 'nrow'
  pairwise : lincon_to_pairwise: warning in matrix(rep(NA,
    (length(t1$fnames) - 1)^2), nc = length(t1$fnames) - 1, nr =
    length(t1$fnames) - 1): correspondance partielle d'argument de 'nc'
    par rapport à 'ncol'
  pairwise : bootstrap: no visible binding for global variable
    'pairwise.t.test'
  pairwise : bootstrap: no visible binding for global variable 'quantile'
  pairwise: no visible global function definition for 'pairwise.t.test'
  pairwise: no visible binding for global variable 'pairwise.t.test'
  pairwise : p.w.t: no visible global function definition for
    'pairwise.wilcox.test'
  pairwise: no visible binding for global variable 'median'
  pairwise : BM_to_pairwise: no visible global function definition for
    'p.adjust'
  pairwise: no visible global function definition for 'median'
  pairwise: warning in matrix(rep(NA, (length(unique_g) - 1)^2), nc =
    (length(unique_g) - 1), nr = (length(unique_g) - 1)): correspondance
    partielle d'argument de 'nr' par rapport à 'nrow'
  pairwise: warning in matrix(rep(NA, (length(unique_g) - 1)^2), nc =
    (length(unique_g) - 1), nr = (length(unique_g) - 1)): correspondance
    partielle d'argument de 'nc' par rapport à 'ncol'
  pairwise: no visible binding for global variable 'quantile'
  pairwise.boot: no visible binding for global variable 'median'
  pairwise.boot: warning in matrix(rep(NA, (length(unique_g) - 1)^2), nc
    = (length(unique_g) - 1), nr = (length(unique_g) - 1)):
    correspondance partielle d'argument de 'nr' par rapport à 'nrow'
  pairwise.boot: warning in matrix(rep(NA, (length(unique_g) - 1)^2), nc
    = (length(unique_g) - 1), nr = (length(unique_g) - 1)):
    correspondance partielle d'argument de 'nc' par rapport à 'ncol'
  pairwise.boot: no visible global function definition for 'median'
  parco: no visible global function definition for 'formula'
  parco: no visible global function definition for 'rainbow'
  pareto: no visible global function definition for 'barplot'
  pareto: no visible global function definition for 'abline'
  pareto: no visible global function definition for 'box'
  pareto: no visible global function definition for 'axis'
  pareto: no visible global function definition for 'mtext'
  pde: no visible global function definition for 'colorRampPalette'
  range_to_sd: no visible global function definition for 'qnorm'
  skewness.norm.test: no visible global function definition for 'rnorm'
  valreg : bptest_lmer: no visible global function definition for
    'residuals'
  valreg : cooks_distance_lmer : <anonymous>: no visible global function
    definition for 'residuals'
  valreg : cooks_distance_lmer : <anonymous>: no visible global function
    definition for 'sigma'
  valreg : count_terms_excluding_intercept: no visible global function
    definition for 'terms'
  valreg: no visible global function definition for 'coef'
  valreg: no visible global function definition for 'formula'
  valreg: no visible binding for global variable 'code'
  valreg: no visible global function definition for 'getCall'
  wilcoxon.cut.test: no visible global function definition for 'median'
  wilcoxon.cut.test: no visible binding for global variable 'k'
  wilcoxon.cut.test: no visible global function definition for 'quantile'
  wilcoxon.cut.test: no visible binding for global variable 'i'
  Undefined global functions or variables:
    BIC IQR abline ansari.test as.dist axis barplot binom.test binomial
    box boxplot chisq.test code coef colorRampPalette complete.cases cov
    cutree delete.response density dfbetas drop1 error exit form formula
    friedman.test games_howell_test getCall get_all_vars glm grid
    hatvalues hclust i influence k logit mahalanobis manova median
    model.response mood.test msg mtext na.pass p.adjust pairwise.t.test
    pairwise.wilcox.test pvals qbinom qchisq qnorm qt quantile rainbow
    rbeta reformulate reg reshape residuals rnorm rstandard rstudent
    setNames sigma str_split_fixed stripchart terms update var verbose
    weight
  Consider adding
    importFrom("grDevices", "colorRampPalette", "rainbow")
    importFrom("graphics", "abline", "axis", "barplot", "box", "boxplot",
               "grid", "mtext", "stripchart")
    importFrom("stats", "BIC", "IQR", "ansari.test", "as.dist",
               "binom.test", "binomial", "chisq.test", "coef",
               "complete.cases", "cov", "cutree", "delete.response",
               "density", "dfbetas", "drop1", "formula", "friedman.test",
               "getCall", "get_all_vars", "glm", "hatvalues", "hclust",
               "influence", "mahalanobis", "manova", "median",
               "model.response", "mood.test", "na.pass", "p.adjust",
               "pairwise.t.test", "pairwise.wilcox.test", "qbinom",
               "qchisq", "qnorm", "qt", "quantile", "rbeta", "reformulate",
               "reshape", "residuals", "rnorm", "rstandard", "rstudent",
               "setNames", "sigma", "terms", "update", "var")
  to your NAMESPACE file.

❯ checking Rd files ... NOTE
  checkRd: (-1) pairwise.Rd:29: Lost braces; missing escapes or markup?
      29 | \item{type}{'mean' for pairwise.t.test(p.adjust.method="holm"), 'median' for pairwise.wilcox.test(p.adjust.method="BH"), 'ks' for ks.test(), 'lincon' for lincon() of {WSR2}}
         |                                                                                                                                                                       ^
  checkRd: (-1) pairwise.Rd:54: Lost braces; missing escapes or markup?
      54 | This function automates the work of the ks.test(), lincon() functions of {WSR2}, pairwise.t.test() and pairwise.wilcox.test() and extracts groups of means or comparisons to a control with the catego() function.
         |                                                                          ^

❯ checking package vignettes ... NOTE
  Package has 'vignettes' subdirectory but apparently no vignettes.
  Perhaps the 'VignetteBuilder' information is missing from the
  DESCRIPTION file?

1 error ✖ | 3 warnings ✖ | 9 notes ✖
Erreur : R CMD check found ERRORs
Exécution arrêtée

Exited with status 1.