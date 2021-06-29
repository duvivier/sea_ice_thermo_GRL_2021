#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')
set vars = ('LHFLX' 'PBLH' 'PSL' 'QFLX' 'QREFHT' 'SHFLX' 'SST' 'TREFHT' 'TS' 'U10' 'PRECC' 'PRECL' 'PRECSC' 'PRECSL' 'PRECT' 'CLDLOW' 'CLDMED' 'CLDHGH' 'FLDS')

##############
# start loops
set c = 1
while ($c <= 2)  # max: 2
set v = 19
while ($v <= 19)  # max: 15

set case = $cases[$c]
set varcode = $vars[$v]

##############
# Input into ncl
##############

        echo 'Creating file of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_spatial_atm_subset.ncl

@ v ++
end
@ c ++
end

