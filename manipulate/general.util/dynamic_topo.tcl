package require topotools

# ==============================================================================                                                                              
# HELP UTILITY COMMAND                                                                                                                                        
# ==============================================================================                                                                              
proc show_help {} {                                                                                                                                           
    puts ""                                                                                                                                                   
    puts "======================================================================"
    puts "VMD DYNAMIC TOPOLOGY UTILITY - HELP MENU"
    puts "======================================================================"
    puts "0. Import topotools"                                                                                                                                
    puts "   package require topotools"                                                                                                                       
    puts ""                                                                                                                                                   
    puts "1. LOAD SCRIPT"                                                                                                                                     
    puts "   source dynamic_topo.tcl"                                                                                                                       
    puts ""                                                                                                                                                   
    puts "2. DISPLAY BONDS on the fly BASED ON DYNAMIC TOPOLOGY MATRIX"
    puts "   usage: dynamic_topology_on frames.top"                     
    puts "======================================================================"
    puts ""
}       

# ------------------------------------------------------------------
# Global state
# ------------------------------------------------------------------
array set ::WBONDS {}        ;# frame -> list of {i j} water bonds (0-based)
set ::SOLUTE_BONDS {}        ;# static bonds not involving water
set ::WMOLID -1
set ::INDEX_BASE 1           ;# set to 1 if the file is Fortran 1-based, 0 if already 0-based

# Automatically display help instructions upon sourcing the file                                                                                              
show_help                                                                                                                                                     

# ------------------------------------------------------------------
# 1. Parse the topology file
# ------------------------------------------------------------------
proc read_water_topology {filename} {
    array unset ::WBONDS
    set fh [open $filename r]
    set frame -1
    set nframes 0

    while {[gets $fh line] >= 0} {
        set line [string trim $line]
        if {$line eq ""} { continue }

        if {[regexp -nocase {^FRAME\s+(\d+)} $line -> f]} {
            set frame $f
            set ::WBONDS($frame) {}
            incr nframes
            continue
        }

        if {$frame < 0} { continue }   ;# data before first FRAME header

        if {[scan $line "%d %d" a b] == 2} {
            # convert to VMD 0-based indices
            lappend ::WBONDS($frame) \
                [list [expr {$a - $::INDEX_BASE}] [expr {$b - $::INDEX_BASE}]]
        }
    }
    close $fh
    puts "read $nframes frames of water connectivity from $filename"
}

# ------------------------------------------------------------------
# 2. Capture the solute bonds once (everything not O–H water)
# ------------------------------------------------------------------
proc capture_solute_bonds {molid water_seltext} {
    set wsel [atomselect $molid $water_seltext]
    set water_idx [$wsel get index]
    $wsel delete

    # fast membership lookup
    array unset iswater
    foreach i $water_idx { set iswater($i) 1 }

    set ::SOLUTE_BONDS {}
    foreach bond [topo getbondlist none -molid $molid] {
        lassign $bond i j
        # keep the bond only if NEITHER atom is a water atom
        if {![info exists iswater($i)] && ![info exists iswater($j)]} {
            lappend ::SOLUTE_BONDS $bond
        }
    }
    puts "kept [llength $::SOLUTE_BONDS] solute bonds (static)"
}

# ------------------------------------------------------------------
# 3. Apply connectivity for one frame
# ------------------------------------------------------------------
proc apply_frame_bonds {frame} {
    if {$::WMOLID < 0} { return }
    if {![info exists ::WBONDS($frame)]} { return }   ;# no data for this frame
    topo setbondlist [concat $::SOLUTE_BONDS $::WBONDS($frame)] -molid $::WMOLID
}

# ------------------------------------------------------------------
# 4. Hook it to frame changes
# ------------------------------------------------------------------
proc dynamic_topology_on {filename {molid top} {water_seltext "resname HOH SOL TIP3 WAT"} {base 1}} {
    if {$molid eq "top"} { set molid [molinfo top] }
    set ::WMOLID $molid
    set ::INDEX_BASE $base

    read_water_topology $filename
    capture_solute_bonds $molid $water_seltext

    # frame numbering: file FRAME 1 -> VMD frame 0
    proc ::_topo_trace {name1 name2 op} {
        global vmd_frame
        apply_frame_bonds [expr {$vmd_frame($::WMOLID) + 1}]
    }
    trace add variable ::vmd_frame($molid) write ::_topo_trace

    apply_frame_bonds [expr {[molinfo $molid get frame] + 1}]
    puts "dynamic topology enabled for molid $molid"
}

proc dynamic_topology_off {} {
    if {$::WMOLID < 0} { return }
    trace remove variable ::vmd_frame($::WMOLID) write ::_topo_trace
    puts "dynamic topology disabled"
}

proc load_dynamic {topfile {pdbfile ""} {water_seltext "resname HOH"} {base 1}} {
    if {$pdbfile ne ""} {
        set molid [mol new $pdbfile type pdb waitfor all]
    } elseif {[llength [molinfo list]] > 0} {
        set molid [molinfo top]
        puts "using already-loaded molecule $molid"
    } else {
        error "no molecule loaded and no pdb given"
    }
    dynamic_topology_on $topfile $molid $water_seltext $base
    return $molid
}
