# KefiR Package Check Summary
# Generated: 2026-02-02 19:15

## Status: 1 ERROR | 6 WARNINGS | 4 NOTES

### CRITICAL ERROR (Must Fix)
1. **m.test**: Formal argument "code" matched by multiple actual arguments
   - Location: Call to `.posthoc_mixed_model()`
   - Issue: `code` parameter appears twice in function call
   - Action: Remove duplicate `code` argument

### WARNINGS (Should Fix for CRAN)
1. **Non-ASCII characters** in 22 R files:
   - Already fixed: sys_posthoc_ANCOVA.R, sys_posthoc_MANOVA.R, sys_normality.R, 
     sys_detect_model_type.R, sys_filter_small_groups.R, sys_formulator_safe.R,
     sys_plot_with_letters.R, sys_mixed_model_analysis.R, sys_diagnostic_influence.R,
     sys_manova_analysis.R (partial)
   - Still need fixing: cooks.distance_lmer.R, evolreg.r, gage_rr.R, m.test.R,
     pairwise.R, pairwise.boot.R, pde.R, sample_bootstrap.R, sys_analyze_ancova_structure.R,
     sys_ancova_analysis.R, sys_auto_preprocess_g.R, sys_control_independence.R,
     sys_formulator.R, sys_multi_factor_analysis.R, sys_one_factor_analysis.R,
     sys_pairwise_medpb2.R, sys_posthoc.R, sys_posthoc_mixed_model.R,
     sys_select_ss_type.R, sys_variance.R, valreg.R

2. **Dependencies in R code**:
   - Missing imports: 'broom', 'dae', 'rcompanion', 'rlang'
   - Missing exports: 'emmeans::cld', 'lmerTest::anova'
   - Unexported objects: 'lmerTest:::as_lmerModLmerTest', 'lmerTest:::summary.lmerModLmerTest'

### NOTES (Good to Fix)
1. **DESCRIPTION meta-information**:
   - Duplicate packages in Depends/Imports/Suggests: 'WRS2', 'tseries'
   
2. **Top-level files** (non-standard):
   - 'INSTALLATION_GUIDE.md', 'agent.md', 'res', 'test_check.bat'
   - Action: Move to .Rbuildignore or remove

3. **Missing documentation**:
   - 'cooks.distance_lmer' needs Rd file
   
4. **Undocumented arguments** in multiple Rd files:
   - catego.Rd: 'debug'
   - control_independence.Rd: 'verbose', 'k'
   - dot-code_ancova.Rd: 'step_num', 'title', 'code_lines'
   - dot-code_anova.Rd: 'step_num', 'title', 'code_lines'
   - dot-formulator.Rd: 'verbose', 'return'
   - dot-manova_analysis.Rd: 'within', 'between'
   - dot-posthoc_MANOVA.Rd: 'code'

5. **Missing NAMESPACE imports** (long list of stats/graphics functions)

## Priority Actions
1. Fix ERROR in m.test (duplicate code argument)
2. Continue ASCII conversion for remaining files
3. Add missing imports to NAMESPACE
4. Document missing arguments
5. Clean up DESCRIPTION duplicates
6. Add .Rbuildignore for non-standard files
