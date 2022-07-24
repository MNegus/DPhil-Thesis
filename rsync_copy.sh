rsync -zarv  --include="*/" --include="*.csv" --include="*.dat" --include="output.txt" --include="log" --exclude="*" polaris:/scratch/negus/stationary_plate_maxlevel_validation . 
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" vanisher:/scratch/negus/imposed_0.00125_maxlevel_validation .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" nocturne:/scratch/negus/imposed_membrane_maxlevel_validation .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" wolverine:/scratch/negus/stationary_plate_maxlevel_validation . 
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" mystique:/scratch/negus/imposed_plate_0.00125_maxlevel_validation .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" shadowcat:/scratch/negus/imposed_membrane_0.05_maxlevel_validation .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" shadowcat:/scratch/negus/imposed_membrane_0.0025_maxlevel_validation .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" mystique:/scratch/negus/updated_imposed_membrane_0.0025_maxlevel_validation .


# Membrane outputs from imposed
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --include="membrane_arr_*.txt" --exclude="*" shadowcat:/scratch/negus/imposed_membrane_0.05_maxlevel_validation/imposed_coeff_0 .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --include="membrane_arr_*.txt" --exclude="*" shadowcat:/scratch/negus/imposed_membrane_0.0025_maxlevel_validation .

# Boundary outputs from imposed
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --include="boundary_output_*.txt" --exclude="*" shadowcat:/scratch/negus/imposed_membrane_0.05_maxlevel_validation/imposed_coeff_0 .
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --include="boundary_output_*.txt" --exclude="*" shadowcat:/scratch/negus/imposed_membrane_0.0025_maxlevel_validation .


# DOUBLE WIDTH
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="log" --exclude="*" micromax:/scratch/negus/DOUBLE_WIDTH_TEST .

# QUADRATIC SUBSTRATE
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="boundary_output_*.txt" --include="log" --exclude="*" micromax:/scratch/negus/imposed_0.00125_quadratic_maxlevel_validation .

# EASY CASE
# rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include="boundary_output_*.txt" --include="log" --exclude="*" nocturne:/scratch/negus/EasyCase_maxlevel_validation .

# rsync -zarv  --include="*/" --include="*.c" --include="*.h" --include="*.csv" --include="output.txt" --include="boundary_output_*.txt" --include="log" --exclude="*" nocturne:/scratch/negus/imposed_0.05_omega_6 .

# IMPOSED 0.05 OMEGA 4
rsync -zarv  --include="*/" --include="*.csv" --include="output.txt" --include='*boundary_outputs/**' --include="*.dat" --include="log" --exclude="*" nocturne:/scratch/negus/imposed_0.05_omega_4_maxlevel_validation .

