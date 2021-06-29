#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')

##############
# start loops
set c = 2
while ($c <= 2)  # max: 2

set case = $cases[$c]

##############
# Input into ncl
##############

        echo 'Creating files for '$case
        ncl 'case          = "'$case'"'\
            ./physics_spatial_atm_turbflux_subset.ncl

        ncl 'case          = "'$case'"'\
            ./physics_spatial_ice_massgain_subset.ncl

        ncl 'case          = "'$case'"'\
            ./physics_spatial_ice_winds_subset.ncl

        ncl 'case          = "'$case'"'\
            ./physics_spatial_ocn_Chl_subset.ncl
@ c ++
end

