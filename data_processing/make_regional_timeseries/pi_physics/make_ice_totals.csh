#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')
set vars = ('aice' 'hi' 'frazil' 'congel' 'snoice' 'dvidtd' 'meltt' 'meltb' 'meltl' 'frazil_plus_congel' 'frazil_plus_congel_plus_snoice' 'totChl_150m' 'totChl_100m' 'totChl_50m')

##############
# start loops
set c = 2
while ($c <= 2)  # max: 2
set v = 1
while ($v <= 14)  # max: 14

set case = $cases[$c]
set varcode = $vars[$v]

##############
# Input into ncl
##############

        echo 'Creating file of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_ice_totals.ncl

@ v ++
end
@ c ++
end

