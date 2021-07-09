#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('TEMP' 'SALT' 'PD')

##############
# start loops
set v1 = 3
while ($v1 <= 3)  # max: 8

set varcode1 = $vars[$v1]

##############
# Input into ncl
##############

        echo 'Creating figures for '$varcode1
        ncl 'varcode1       = "'$varcode1'"'\
            ./plot_transect_spatial_mean_diff_hmxl_winter.ncl
        
        rm *.ps

@ v1 ++
end


