#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')
set vars = ('Tinz' 'Sinz')

##############
# start loops
set c = 1
while ($c <= 2)  # max: 2
set v = 1
while ($v <= 2)  # max: 3

set case = $cases[$c]
set varcode = $vars[$v]

##############
# Input into ncl
##############

        echo 'Creating file of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_avg_to_2d_ice_prop_subset.ncl

#physics_spatial_vert_ice_prop_subset.ncl

@ v ++
end
@ c ++
end


