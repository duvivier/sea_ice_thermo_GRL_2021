#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set cases = ('b.e21.B1850.f09_g17.CMIP6-piControl.001b' 'b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b')
set vars = ('LHFLX' 'PBLH' 'PSL' 'QFLX' 'QREFHT' 'SHFLX' 'SST' 'TREFHT' 'TS' 'U10' 'PRECC' 'PRECL' 'PRECSC' 'PRECSL' 'PRECT' 'TURBFLX')

##############
# start loops
set c = 2
while ($c <= 2)  # max: 2
set v = 1
while ($v <= 16)  # max: 16

set case = $cases[$c]
set varcode = $vars[$v]

##############
# Input into ncl
##############

        echo 'Making regional averages of '$varcode
        ncl 'varcode       = "'$varcode'"'\
            'case          = "'$case'"'\
            ./physics_atm_avgs.ncl

@ v ++
end
@ c ++
end

