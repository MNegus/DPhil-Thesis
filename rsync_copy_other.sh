# Stationary plate data
#rsync -zarv polaris:/scratch/negus/stationary_plate_maxlevel_validation/axi_1/max_level_13/ Stationary_Plate/axi

# Moving plate test data
#rsync -zarv nocturne:/scratch/negus/ParamTest_3_8_2022/ALPHA-2_BETA-0_GAMMA-20 Moving_Plate/

# Plate parameter sweep (just keep log and turnover files)
#rsync -zarv  --include="*/" --include="*.csv" --include="*.dat" --include="output.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" micromax:/scratch/negus/Plate_Parameter_Runs .

# Plate parameter sweep (just keep log and turnover files)
#rsync -zarv  --include="*/" --include="*.csv" --include="*.dat" --include="output.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" micromax:/scratch/negus/LogSpacedParameterRuns .


# Stationary and moving plate test data with increased output
#rsync -zarv --include="*/" --include="*.csv" --include="*.dat" --include="output.txt" --include="interface_*.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" polaris:/scratch/negus/stationary_plate_output/ Stationary_Plate_Output/axi
#rsync -zarv --include="*/" --include="*.csv" --include="*.dat" --include="output.txt" --include="interface_*.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" polaris:/scratch/negus/moving_plate_output/ Moving_Plate_Output

# Rubber run with membrane arrays
#rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat" --include="membrane_arr_*.txt" --include="interface_*.txt" --include="log" --include="*.c" --include="*.h" --exclude="*"   nightcrawler:/scratch/negus/RubberRuns/ALPHA_1.1_GAMMA_668 RubberRuns/

# Stationary membrane case (until the better one finishes)
#rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat" --include="interface_*.txt" --include="membrane_arr_*.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" nocturne:/scratch/negus/MembraneStationaryCase/* Stationary_Membrane 

# Membrane parameter runs
rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat" --include="membrane_arr_2250.txt" --include="membrane_arr_3250.txt" --include="membrane_arr_7000.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" micromax:/scratch/negus/MembraneParameterRuns/GAMMA_varying MembraneParameterRuns/

rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat" --include="membrane_arr_2250.txt" --include="membrane_arr_3250.txt" --include="membrane_arr_7000.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" wolverine:/scratch/negus/MembraneParameterRuns/DELTA_varying MembraneParameterRuns/

rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat"  --include="membrane_arr_2250.txt" --include="membrane_arr_3250.txt" --include="membrane_arr_7000.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" nocturne:/scratch/negus/MembraneParameterRuns/BETA_varying MembraneParameterRuns/

rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat"  --include="membrane_arr_2250.txt" --include="membrane_arr_3250.txt" --include="membrane_arr_7000.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" micromax:/scratch/negus/MembraneParameterRuns/NEW_BETA_varying MembraneParameterRuns/

rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat"  --include="membrane_arr_2250.txt" --include="membrane_arr_3250.txt" --include="membrane_arr_7000.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" micromax:/scratch/negus/MembraneParameterRuns/NEW_GAMMA_varying MembraneParameterRuns/

# DNS validation
rsync -zarc --include="*/" --include="*.csv" --include="output.txt" --include="*.dat"  --include="membrane_arr_2250.txt" --include="membrane_arr_3250.txt" --include="membrane_arr_7000.txt" --include="log" --include="*.c" --include="*.h" --exclude="*" swamp-beast:/scratch/negus/DNS_Chapter_validation/imposed_0.05_omega_4_maxlevel_validation .



