package require topotools

# ==============================================================================
# HELP UTILITY COMMAND
# ==============================================================================
proc show_help {} {
    puts ""
    puts "================================================================================"
    puts "VMD TOPOTOOLS BOND UTILITY - HELP MENU"
    puts "================================================================================"
    puts "0. Import topotools"
    puts "   package require topotools"
    puts ""
    puts "1. LOAD SCRIPT"
    puts "   source vmd_bond_utils.tcl"
    puts ""
    puts "2. SCAN AND IDENTIFY LONG BONDS"
    puts "   usage: find_long_bonds <AtomType1> <AtomType2> <Cutoff_Angstroms>"
    puts ""
    puts "3. DELETE LONG BONDS"
    puts "   usage: delete_long_bonds <AtomType1> <AtomType2> <Cutoff_Angstroms>"
    puts ""
    puts "4. FIND SHORT NONBONDS"
    puts "   usage: find_short_nonbonds <AtomType1> <AtomType2> <Cutoff_Angstroms>"
    puts "================================================================================"
    puts ""
}

# ================================================================================
#                            TIP
#    If your topology relies on unique atom names instead of types,
#    change 'get type' to 'get name' inside this script file.
# ================================================================================

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


# ==============================================================================
# FUNCTION 3: FIND SHORT NONBONDS
# ==============================================================================
proc find_short_nonbonds {type1 type2 cutoff} {
    # 1. Select all atoms matching the requested types
    set sel_all1 [atomselect top "type $type1"]
    set sel_all2 [atomselect top "type $type2"]
    
    set count 0
    puts "\n--- Searching for nonbonded $type1-$type2 pairs closer than $cutoff Å ---"
    
    # 2. Use VMD's optimized 'measure contacts' to find pairs within the cutoff radius
    # This returns a list of two sub-lists: { {indices from sel1} {indices from sel2} }
    set contacts [measure contacts $cutoff $sel_all1 $sel_all2]
    set list1 [lindex $contacts 0]
    set list2 [lindex $contacts 1]
    
    # Clean up the initial heavy selections
    $sel_all1 delete
    $sel_all2 delete
    
    # 3. Get the existing covalent bond matrix to filter out bonded pairs
    # 'topo getbondlist' provides pairs like { {0 1} {4 5} ... }
    set bond_list [topo getbondlist none]
    
    # Convert bond list into an array/dict for instant O(1) lookups
    foreach bond $bond_list {
        set b1 [lindex $bond 0]
        set b2 [lindex $bond 1]
        set bonded_matrix($b1,$b2) 1
        set bonded_matrix($b2,$b1) 1
    }
    
    # 4. Loop through the spatial contacts and isolate nonbonded pairs
    foreach id1 $list1 id2 $list2 {
        # Skip the pair if they are the exact same atom (relevant if type1 == type2)
        if {$id1 == $id2} { continue }
        
        # Skip the pair if a chemical bond already connects them
        if {[info exists bonded_matrix($id1,$id2)]} { continue }
        
        # Prevent double-counting the same pair in different orders (relevant if type1 == type2)
        if {$type1 == $type2 && $id1 > $id2} { continue }
        
        # Calculate the exact distance
        set s1 [atomselect top "index $id1"]
        set s2 [atomselect top "index $id2"]
        
        set coord1 [lindex [$s1 get {x y z}] 0]
        set coord2 [lindex [$s2 get {x y z}] 0]
        set dist [veclength [vecsub $coord1 $coord2]]
        
        set name1 [$s1 get name]
        set name2 [$s2 get name]
        
        puts "Index $id1 ($name1) -- Index $id2 ($name2) | Distance: [format "%.2f" $dist] Å"
        incr count
        
        $s1 delete
        $s2 delete
    }
    
    puts "--------------------------------------------------------"
    puts "Total short nonbonded $type1-$type2 pairs found: $count\n"
}
