#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')
set vars = ('frazil' 'congel' 'snoice' 'frazil_plus_congel' 'frazil_plus_congel_plus_snoice' 'meltt' 'meltb' 'meltl' 'daidtd' 'daidtt' 'dvidtd' 'dvidtt' 'fswdn' 'fswthru' 'aice' 'hi' 'uatm' 'vatm' 'wsatm' 'fresh' 'fsalt' 'fhocn' 'HMXL' 'HMXL_DR' 'HMXL_DR2' 'XMXL' 'XMXL_DR' 'PD_50m' 'PD_100m' 'PD_200m' 'PD_300m' 'PD_500m' 'TEMP_50m' 'TEMP_100m' 'TEMP_200m' 'TEMP_300m' 'TEMP_500m' 'SALT_50m' 'SALT_100m' 'SALT_200m' 'SALT_300m' 'SALT_500m' 'totChl')

##############
# start loops
set c = 1
while ($c <= 2)  # max: 2
set v = 43
while ($v <= 43)  # max: 43

set case = $cases[$c]
set varcode = $vars[$v]

##############
# Input into ncl
##############

        echo 'Making seasonal averages of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_spatial_ice_ocn_seasonal_avgs.ncl

        echo 'Making monthly averages of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_spatial_ice_ocn_monthly_avgs.ncl

@ v ++
end
@ c ++
end

