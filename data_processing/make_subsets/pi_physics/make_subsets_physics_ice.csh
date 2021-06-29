#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')

set vars = ('aice' 'daidtd' 'daidtt' 'dvidtd' 'dvidtt' 'fswdn' 'fswthru' 'hi' 'ice_present' 'frazil' 'congel' 'snoice' 'meltt' 'meltb' 'meltl' 'uatm' 'vatm' 'fresh' 'fsalt' 'fhocn')
#set vars = ('aice' 'daidtd' 'daidtt' 'dvidtd' 'dvidtt' 'fswdn' 'fswthru' 'hi' 'ice_present' 'sidmassgrowthwat' 'sidmassgrowthbot' 'sidmasssi' 'sidmassdyn' 'sidmassmelttop' 'sidmassmeltbot' 'sidmasslat' 'sidmassevapsubl' 'uatm' 'vatm')

##############
# start loops
set c = 2
while ($c <= 2)  # max: 2
set v = 1
while ($v <= 20)  # max: 20

set case = $cases[$c]
set varcode = $vars[$v]

##############
# Input into ncl
##############

        echo 'Creating file of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_spatial_ice_subset.ncl

@ v ++
end
@ c ++
end


