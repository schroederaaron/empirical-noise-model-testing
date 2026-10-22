# power_constants.R
# Constants shared by power_test.R and estimate_dispersion_trend.R (note D WP2).
# POWER_LIB_SIZE is the simulator's library size AND the depth at which the P+
# dispersion trend's mu is defined; a parameter file fitted at another depth is
# rejected by power_test.R (check_disp_params).
POWER_LIB_SIZE <- 4e7
