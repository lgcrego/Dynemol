###################################################################
#Script Name    : MO-show.tcl                                                                                             
#Description    : vmd -e MO-show.tcl                                                                               
#Args           : no arguments                                                                                          
#Author         : Luis G C Rego
#date           : 12/Oct/2023
###################################################################

# Display settings
display projection   Orthographic
display nearclip set 0.000000
display farclip  set 10.000000
display depthcue   off
# Set the directory where your cube files are located
set directory_path "./"

# Get a list of all cube files in the directory and sort them
set file_list [lsort [glob -directory $directory_path -nocomplain *.cube]]

# Use the first file in the list to create a molecule
set updmol [mol new [lindex $file_list 0] type cube waitfor all]

# Loop through the rest of the files and use mol addfile
for {set i 1} {$i < [llength $file_list]} {incr i} {
    mol addfile [lindex $file_list $i] type cube waitfor all
}

mol delrep 0 top
mol representation CPK 0.500000 0.300000 20.000000 16.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep top
mol color Name
mol addrep top
mol representation Isosurface -0.03 0.0 0.0 0.0
mol color ColorID 1
mol selection {all}
mol addrep top
mol representation Isosurface 0.03 0.0 0.0 0.0
mol color ColorID 0
mol selection {all}

# store name of the isosurface representation (id=3) for later use
mol addrep top
set updrep_2 [mol repname top 2]
set updrep_3 [mol repname top 3]
mol rename top {Golden Shot}

# use the volumetric data set for the isosurface corresponding to the frame.
# $updmol contains the id of the molecule and $updrep the (unique) name of 
# the isosurface representation
#
proc update_iso {args} { 
    global updmol
    global updrep_2
    global updrep_3

    # get representation id and return if invalid
    set repid [mol repindex $updmol $updrep_3]
    if {$repid < 0} { return }

    # update representation but replace the data set 
    # id with the current frame number.
    set frame [molinfo $updmol get frame]
    lassign [molinfo $updmol get "{rep $repid}"] rep
    mol representation [lreplace $rep 2 2 $frame]
    mol color ColorID 0
    mol modrep $repid $updmol

    # get representation id and return if invalid
    set repid [mol repindex $updmol $updrep_2]
    if {$repid < 0} { return }

    # update representation but replace the data set 
    # id with the current frame number.
    set frame [molinfo $updmol get frame]
    lassign [molinfo $updmol get "{rep $repid}"] rep
    mol representation [lreplace $rep 2 2 $frame]
    mol color ColorID 1
    mol modrep $repid $updmol

}

trace variable vmd_frame($updmol) w update_iso
animate goto 0

