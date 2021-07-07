#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('LHFLX' 'SHFLX')

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
            ./plot_spatial_diffs_winter.ncl
        
        rm *.ps

@ v1 ++
end


