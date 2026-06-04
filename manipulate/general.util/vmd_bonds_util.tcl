package require topotools

# ==============================================================================
# HELP UTILITY COMMAND
# ==============================================================================
proc show_help {} {
    puts ""
    puts "================================================================================"
    puts "VMD TOPOTOOLS BOND UTILITY - HELP MENU"
    puts "================================================================================"
    puts "0. Import topotools; use VMD's TkConsole"
    puts "   package require topotools"
    puts ""
    puts "1. LOAD SCRIPT"
    puts "   source vmd_bond_utils.tcl"
    puts ""
    puts "2. SCAN AND IDENTIFY LONG BONDS (Read-Only)"
    puts "   find_long_bonds <AtomType1> <AtomType2> <Cutoff_Angstroms>"
    puts "   Example: find_long_bonds H O 2.5"
    puts ""
    puts "3. DELETE LONG BONDS (Modifies Topology)"
    puts "   delete_long_bonds <AtomType1> <AtomType2> <Cutoff_Angstroms>"
    puts "   Example: delete_long_bonds C C 1.8"
    puts ""
    puts "4. PRO-TIP"
    puts "   If your topology relies on unique atom names instead of types,"
    puts "   change 'get type' to 'get name' inside this script file."
    puts "================================================================================"
    puts ""
}

# Automatically display help instructions upon sourcing the file
show_help


# ==============================================================================
# FUNCTION 1: IDENTIFY LONG BONDS
# ==============================================================================
proc find_long_bonds {type1 type2 cutoff} {
    set bond_list [topo getbondlist none]
    set count 0
    
    puts "\n--- Searching for $type1-$type2 bonds longer than $cutoff Å ---"
    
    foreach bond $bond_list {
        set id1 [lindex $bond 0]
        set id2 [lindex $bond 1]
        
        set sel1 [atomselect top "index $id1"]
        set sel2 [atomselect top "index $id2"]
        
        set t1 [$sel1 get type]
        set t2 [$sel2 get type]
        
        if { ($t1 == $type1 && $t2 == $type2) || ($t1 == $type2 && $t2 == $type1) } {
            set coord1 [lindex [$sel1 get {x y z}] 0]
            set coord2 [lindex [$sel2 get {x y z}] 0]
            set dist [veclength [vecsub $coord1 $coord2]]
            
            if { $dist > $cutoff } {
                set name1 [$sel1 get name]
                set name2 [$sel2 get name]
                puts "Index $id1 ($name1) -- Index $id2 ($name2) | Distance: [format "%.2f" $dist] Å"
                incr count
            }
        }
        $sel1 delete
        $sel2 delete
    }
    puts "--------------------------------------------------------"
    puts "Total overlong $type1-$type2 bonds found: $count\n"
}


# ==============================================================================
# FUNCTION 2: DELETE LONG BONDS
# ==============================================================================
proc delete_long_bonds {type1 type2 cutoff} {
    set bond_list [topo getbondlist none]
    set bonds_to_delete {}
    
    puts "\n--- Scanning for $type1-$type2 bonds longer than $cutoff Å ---"
    
    foreach bond $bond_list {
        set id1 [lindex $bond 0]
        set id2 [lindex $bond 1]
        
        set sel1 [atomselect top "index $id1"]
        set sel2 [atomselect top "index $id2"]
        
        set t1 [$sel1 get type]
        set t2 [$sel2 get type]
        
        if { ($t1 == $type1 && $t2 == $type2) || ($t1 == $type2 && $t2 == $type1) } {
            set coord1 [lindex [$sel1 get {x y z}] 0]
            set coord2 [lindex [$sel2 get {x y z}] 0]
            set dist [veclength [vecsub $coord1 $coord2]]
            
            if { $dist > $cutoff } {
                set name1 [$sel1 get name]
                set name2 [$sel2 get name]
                puts "Targeted for deletion: Index $id1 ($name1) -- Index $id2 ($name2) | Distance: [format "%.2f" $dist] Å"
                lappend bonds_to_delete [list $id1 $id2]
            }
        }
        $sel1 delete
        $sel2 delete
    }
    
    set delete_count [llength $bonds_to_delete]
    if { $delete_count > 0 } {
        foreach bond $bonds_to_delete {
            set b1 [lindex $bond 0]
            set b2 [lindex $bond 1]
            ::TopoTools::delbond top $b1 $b2
        }
        mol reanalyze top
        puts "--------------------------------------------------------"
        puts "Successfully deleted $delete_count overlong $type1-$type2 bonds."
    } else {
        puts "--------------------------------------------------------"
        puts "No overlong $type1-$type2 bonds found. Nothing to delete."
    }
    puts ""
}

