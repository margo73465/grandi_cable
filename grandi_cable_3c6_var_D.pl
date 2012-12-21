#!/usr/bin/perl

#This is a batch file that uses the Perl language
#For Perl syntax, see the book "Programming Perl"
#
#To run this, just type perl grandi_cable_3c6_var_D.pl
#
#Note that the first line of this file must remain unchanged (#!/usr/bin/perl) so
#that the system knows that this file should be processed using the perl interpreter
#(thus, the file extension does not need to be .pl ... I use that as my standard)


$screen_log_file = "grandi_cable_3c6_varD.log";
unlink $screen_log_file;


system "gcc grandi_cable_3c6_varD.c -lm -O3 -fomit-frame-pointer -fopenmp -D __ompEnabled -o go1 >> $screen_log_file";
system "echo >> $screen_log_file";

$D = 0.0014;
	
$outputI = "grandiI1_cable_3c6_1Hz_def_D1c4_dx0c0125_dtmax0c01.dat"; 
unlink $OutputI;
print "outputI = $outputI\n";

$outputI2 = "grandiI2_cable_3c6_1Hz_def_D1c4_dx0c0125_dtmax0c01.dat";
unlink $OutputI;
print "outputI2 = $outputI2\n";
	
$outputAPD = "grandi_APD_cable_3c6_1Hz_def_D1c4_dx0c0125_dtmax0c01.dat";  
unlink $outputAPD;
print "outputAPD = $outputAPD\n";
	
$outputVdmax = "grandi_Vdmax_cable_3c6_1Hz_def_D1c4_dx0c0125_dtmax0c01.dat"; 
unlink $outputVdmax;
print "outputVdmax = $outputVdmax\n";

$outputCV = "grandi_CV_cable_3c6_1Hz_def_D1c4_dx0c0125_dtmax0c01.dat";
unlink $outputCV;
print "outputCV = $outputCV\n";

system "./go1 $D $outputI $outputI2 $outputAPD $outputVdmax $outputCV >> $screen_log_file";
system "echo >> $screen_log_file";

print "Simulation with D = $D completed.\n"; 



$D = 0.0004;
	
$outputI = "grandiI1_cable_3c6_1Hz_def_D0c4_dx0c0125_dtmax0c01.dat"; 
unlink $OutputI;
print "OutputI = $OutputI\n";

$outputI2 = "grandiI2_cable_3c6_1Hz_def_D0c4_dx0c0125_dtmax0c01.dat";
unlink $OutputI;
print "OutputI2 = $OutputI2\n";
	
$outputAPD = "grandi_APD_cable_3c6_1Hz_def_D0c4_dx0c0125_dtmax0c01.dat";  
unlink $outputAPD;
print "outputAPD = $outputAPD\n";
	
$outputVdmax = "grandi_Vdmax_cable_3c6_1Hz_def_D0c4_dx0c0125_dtmax0c01.dat"; 
unlink $outputVdmax;
print "outputVdmax = $outputVdmax\n";

$outputCV = "grandi_CV_cable_3c6_1Hz_def_D0c4_dx0c0125_dtmax0c01.dat";
unlink $outputCV;
print "outputCV = $outputCV\n";

system "./go1 $D $outputI $outputI2 $outputAPD $outputVdmax $outputCV >> $screen_log_file";
system "echo >> $screen_log_file";

print "Simulation with D = $D completed.\n"; 



$D = 0.00015;
	
$outputI = "grandiI1_cable_3c6_1Hz_def_D0c15_dx0c0125_dtmax0c01.dat"; 
unlink $OutputI;
print "OutputI = $OutputI\n";

$outputI2 = "grandiI2_cable_3c6_1Hz_def_D0c15_dx0c0125_dtmax0c01.dat";
unlink $OutputI;
print "OutputI2 = $OutputI2\n";
	
$outputAPD = "grandi_APD_cable_3c6_1Hz_def_D0c15_dx0c0125_dtmax0c01.dat";  
unlink $outputAPD;
print "outputAPD = $outputAPD\n";
	
$outputVdmax = "grandi_Vdmax_cable_3c6_1Hz_def_D0c15_dx0c0125_dtmax0c01.dat"; 
unlink $outputVdmax;
print "outputVdmax = $outputVdmax\n";

$outputCV = "grandi_CV_cable_3c6_1Hz_def_D0c15_dx0c0125_dtmax0c01.dat";
unlink $outputCV;
print "outputCV = $outputCV\n";

system "./go1 $D $outputI $outputI2 $outputAPD $outputVdmax $outputCV >> $screen_log_file";
system "echo >> $screen_log_file";

print "Simulation with D = $D completed.\n"; 
