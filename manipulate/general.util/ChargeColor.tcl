set mol [mol new frames.pdb waitfor all] 
set sel [atomselect $mol all] 
set nf [molinfo $mol get numframes] 
set atm [molinfo $mol get numatoms]
set fp [open NetCharge.inpt r] 
set line ""
for {set i 0} {$i < $nf} {incr i} { 
  gets $fp line 
  $sel frame $i 
  $sel set user $line 
  set fixed_min_value -1.0
  set fixed_max_value 1.0
  mol scaleminmax 0 0 $fixed_min_value $fixed_max_value
} 
close $fp 
$sel delete 
display projection orthographic
display depthcue off
color scaleoffset 0.08
animate goto 0
mol modcolor 0 0 User
