#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('aice' 'hi' 'congel' 'frazil' 'snoice' 'frazil_plus_congel_plus_snoice' 'PSL' 'PBLH' 'U10' 'PRECT' 'TREFHT' 'QREFHT' 'TURBFLX' 'LHFLX' 'SHFLX' 'TEMP_100m' 'SALT_100m' 'PD_100m' 'HMXL' 'fresh' 'totChl_150m' 'totChl_100m' 'totChl_50m')

##############
# start loops
set v1 = 1
while ($v1 <= 2)  # max: 23

set varcode1 = $vars[$v1]

##############
# Input into ncl
##############

        echo 'Creating figures for '$varcode1
        
        ncl 'varcode1       = "'$varcode1'"'\
            ./plot_spatial_diffs_winter_unmasked.ncl

        ncl 'varcode1       = "'$varcode1'"'\
            ./plot_spatial_diffs_summer_unmasked.ncl

        rm *.ps

@ v1 ++
end


