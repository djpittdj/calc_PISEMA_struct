package require bigdcd

proc calc_PISEMA {frame} {
	global num_frames z_axis delta
	global DC_sim_list CS_sim_list resid_lst B0PAF_lst
	global nu0 sigma11 sigma22 sigma33

	foreach res $resid_lst {
		set atomC [atomselect top "protein and resid [expr {$res-1}] and name C"]
		set atomN [atomselect top "protein and resid $res and name N"]
		set atomH [atomselect top "protein and resid $res and name HN"]
		set vecC [measure center $atomC]
		set vecN [measure center $atomN]
		set vecH [measure center $atomH]
		set vecNH [vecnorm [vecsub $vecH $vecN]]
		set vecNC [vecnorm [vecsub $vecC $vecN]]

		# dipolar coupling
		set zNH [lindex $vecNH 2]
		set DC_sim [expr {0.5 * $nu0 * (3.0 * $zNH * $zNH - 1.0)}]
		set DC_sim_list($res) [expr {$DC_sim_list($res) + $DC_sim}]

		# chemical shift
		set e2 [vecnorm [veccross $vecNC $vecNH]]
		set mat [transabout $e2 $delta]
		set vecNH4 [lappend vecNH 0.0]
		set e3 [vecnorm [vectrans $mat $vecNH4]]
		set e3 [lrange $e3 0 2]
		set e1 [vecnorm [veccross $e2 $e3]]
		set e1z [lindex $e1 2]
		set e2z [lindex $e2 2]
		set e3z [lindex $e3 2]

		set CS_sim [expr {$sigma11 * $e1z * $e1z + $sigma22 * $e2z * $e2z + $sigma33 * $e3z * $e3z}]
		set CS_sim_list($res) [expr {$CS_sim_list($res) + $CS_sim}]

		set e1 [lappend e1 0]
		set e2 [lappend e2 0]
		set e3 [lappend e3 0]
		set M [list $e1 $e2 $e3 [list 0 0 0 1]]
		set B0PAF_lst($res) [coordtrans $M $z_axis]

		$atomC delete
		$atomN delete
		$atomH delete
	}
	incr num_frames
}

set num_frames 0
set sigma11 57.3
set sigma22 81.2
set sigma33 227.8
set nu0 10.735
set z_axis [list 0 0 -1.0]
set DEG2RAD [expr {3.14159265358979323846/180.0}]
set delta -17.

set resid_lst [list 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]

array set DC_sim_list {}
foreach res $resid_lst {
	set DC_sim_list($res) 0.
}

array set CS_sim_list {}
foreach res $resid_lst {
	set CS_sim_list($res) 0.
}

array set B0PAF_lst {}
foreach res $resid_lst {
	set B0PAF_lst($res) []
}

mol new ideal_Ala16.tilt30.pdb

bigdcd calc_PISEMA ideal_Ala16.tilt30.pdb
bigdcd_wait

set outStream_dc [open dc_calc.xvg w]
set outStream_cs [open cs_calc.xvg w]
set outStream_pisema [open pisema_calc.xvg w]
set outStream_B0PAF [open B0PAF_struct.xvg w]
foreach res $resid_lst {
	set DC_sim_list($res) [expr {$DC_sim_list($res) / double($num_frames)}]
	set CS_sim_list($res) [expr {$CS_sim_list($res) / double($num_frames)}]
	puts $outStream_dc [format "%8i%30.20f" $res $DC_sim_list($res)]
	puts $outStream_cs [format "%8i%30.20f" $res $CS_sim_list($res)]
	puts $outStream_pisema [format "%8i%30.20f%30.20f" $res $CS_sim_list($res) $DC_sim_list($res)]
	puts -nonewline $outStream_B0PAF [format "%8i" $res]
	set B0PAF $B0PAF_lst($res)
	for {set i 0} {$i < [llength $B0PAF]} {incr i} {
		puts -nonewline $outStream_B0PAF [format "%30.20f" [lindex $B0PAF $i]]
	}
	puts $outStream_B0PAF ""
}
close $outStream_dc
close $outStream_cs
close $outStream_pisema
close $outStream_B0PAF

exit
