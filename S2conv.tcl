# --------------------------------------------------------------------------------------------------
# 2D Cantilever Column -- Build Model
# nonlinearBeamColumn element, uniaxial inelastic section
# Units are N,mm
#
#    ^Y
#    |
#    2       ___
#    |        | 
#    |        |
#    |        |
#  (1)        L
#    |        |
#    |        |
#    |        |
#   =1=      _|_  -------->X
#

# SET UP ----------------------------------------------------------------------------
wipe;					# clear memory of all past model definitions
model BasicBuilder -ndm 2 -ndf 3;	# Define the model builder, ndm=#dimension, ndf=#dofs
set dataDir Data;			# set up name for data directory
#file mkdir $dataDir/; 			# create data directory

set PI [expr 2*asin(1.0)];
set g 9810; 	# gravitational acceleration
	
# define GEOMETRY -------------------------------------------------------------
set L 1600; 			# column length
set H 400;			# Column Depth
set B 400;			# Column Width

# calculated mass and axial load
set Q 2112000; # Imposed axial load
set unitwt 	24.0e-6;
set Weight [expr $unitwt*$L*$H*$B];		# Self-weight
set P [expr -$Weight-$Q];		# Total axial load
set Mass [expr $Weight/$g];		# nodal mass
# calculated geometry parameters
set A [expr $B*$H];			# cross-sectional area
set Iz [expr 1./12.*$B*pow($H,3)]; 	# Column moment of inertia
#puts "weight is $Weight"

set nel 4
set nnode [expr $nel + 1]
set elsize [expr $L/$nel]

# nodal coordinates:
#node 1 0 0;				# node#, X, Y
#node 2 0 $L;

for {set x 0} {$x < $nnode } {incr x} {
   set coord [expr $x*$elsize]
   set num [expr $x + 1]
   node $num 0 $coord
   puts "$num  0 $coord"
} 		

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 				# node DX DY RZ

# we need to set up parameters that are particular to the model.
set IDctrlNode $nnode;			# node where displacement is read for displacement control
set IDctrlDOF 1;			# degree of freedom of displacement read for displacement control
set iSupportNode "1";			# define support node, if needed.

# nodal masses:
mass $nnode $Mass 0. 0.;			# node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes

# Define ELEMENTS & SECTIONS -------------------------------------------------------------
set SecTag 1;					# assign a tag number to the column section	
# define section geometry
set cover 17;    # column cover till CL of transverse reinforcement
set steellayer 29;			# Column cover to reinforcing steel NA.
set numBarsLayer1 4;				# number of longitudinal-reinforcement bars in column, layer 1 
set numBarsLayer2 2;				# number of longitudinal-reinforcement bars in column, layer 2
set numBarsLayer3 2;				# number of longitudinal-reinforcement bars in column, layer 3
set numBarsLayer4 4;				# number of longitudinal-reinforcement bars in column, layer 4

set BarDiam 16;			# diameter of corner bars	
set BarArea [expr $PI*pow($BarDiam/2,2)]; 	# area of one corner bar

# MATERIAL parameters -------------------------------------------------------------------
set IDconcCore 1; 				# material ID tag -- confined core concrete
set IDconcCover 2;
set IDreinf 3;
set IDcMinMax 4;
set IDccMinMax 5;
set IDsMinMax 6;


# CONCRETE -----------------------------------------------------------------------------
# nominal concrete compressive strength
set fc 44;				# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
# unconfined concrete
set fc1 $fc;					# UNCONFINED concrete, maximum stress
set ec1 0.002;			# strain at maximum strength of unconfined concrete
#set Ec [expr 2.0*$fc1/$ec1];		# Concrete Elastic Modulus
set Ec [expr 4700*sqrt($fc)];
puts "$Ec"
set Ed -7333.3333;

# confined concrete
set fcc1 48.32;
set ecc1 0.0021964;				# strain at peak stress 
#set Ecc [expr 2.0*$fc1/$ec1];
set Ecc [expr 4700*sqrt($fc)];
set Edc -1181.9683;


# -----------
#uniaxialMaterial uniaxialJ2Plasticity    tag       E   fc  ec0  Ed              
# Core concrete (confined)
uniaxialMaterial NLConcrete $IDconcCore $Ecc $fcc1 $ecc1 $Edc;	# build core concrete (confined)

# Cover concrete (unconfined)
uniaxialMaterial NLConcrete $IDconcCover $Ec $fc1 $ec1 $Ed;	# build cover concrete (unconfined)



# STEEL
# Reinforcing steel 
set Fy 446;				# STEEL yield stress
set Es 200000;			# modulus of steel
set Bs 0.019; 					# strain-hardening ratio 
set Hs [expr $Bs*$Es];
set R0 15;					# control the transition from elastic to plastic branches
set cR1 0.925;					# control the transition from elastic to plastic branches
set cR2 0.15;					# control the transition from elastic to plastic branches

set ets 0.069598;
set ecs -0.069598;

# -----------   

#uniaxialMaterial UniaxialJ2Plasticity $IDreinf $Es $Fy 0.0 $Hs;

uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2;

#uniaxialMaterial MinMax $IDsMinMax $IDreinf -min $ecs -max $ets
	


puts "All material variables have been defined"
# ------------------------------------------

# FIBER SECTION properties -------------------------------------------------------------
# symmetric section
#                   y
#                   ^
#                   |     
#            --------------------     --   --
#           |  o    o    o  |       -- cover
#                               
#  z <---   |  o    +    o  |     
#     
#           |  o    o    o  |       -- cover
#            --------------------     --   --
#           |-------B-------|
#
#                     y
#                     ^
#                     |    
#          ---------------------
#          |\      cover      /|
#          | \------Top------/ |
#          |c|               |c|
#          |o|               |o|
#  z <-----|v|      core     |v|  H
#          |e|               |e|
#          |r|               |r|
#          | /------Bot------\ |
#          |/      cover      \|
#          ---------------------
#                    B
#
# RC section: 
   set coverY [expr $H/2.0];		# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coverZ [expr $B/2.0];		# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coreY [expr $coverY-$cover];     # The distance from the section z-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
   set coreZ [expr $coverZ-$cover];     # The distance from the section y-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
   set steeldistY [expr $coverY-$steellayer];
   set steeldistZ [expr $coverZ-$steellayer];
   set nfCoreY 12;			# number of fibers for concrete in y-direction - core concrete
   set nfCoreZ 1;			# number of fibers for concrete in z-direction - 
   set nfCoverY 12;			# number of fibers for concrete in y-direction -- cover concrete
   set nfCoverZ 1;			# number of fibers for concrete in z-direction
 	

section NLFiber $SecTag {;		# Define the fiber section
	# Define the core patch 		   $yI    $zI     $yJ     $zJ    $yK     $zK    $yL    $zL
	patch quadr $IDconcCore $nfCoreZ $nfCoreY -$coreY $coreZ -$coreY -$coreZ $coreY -$coreZ $coreY $coreZ;

	# Define the four cover patches       $yI     $zI      $yJ    $zJ    $yK    $zK    $yL     $zL
	patch quadr $IDconcCover 1 $nfCoverY -$coverY $coverZ -$coreY $coreZ $coreY $coreZ $coverY $coverZ;
	patch quadr $IDconcCover 1 $nfCoverY -$coreY -$coreZ -$coverY -$coverZ $coverY -$coverZ $coreY -$coreZ;
	patch quadr $IDconcCover $nfCoverZ 1 -$coverY $coverZ -$coverY -$coverZ -$coreY -$coreZ -$coreY $coreZ;
	patch quadr $IDconcCover $nfCoverZ 1  $coreY $coreZ $coreY -$coreZ $coverY -$coverZ $coverY $coverZ;

	# Define reinforcement layers
	layer straight $IDreinf $numBarsLayer1 $BarArea -$steeldistY $steeldistZ -$steeldistY -$steeldistZ;				# bottom layer reinfocement
	#layer straight $IDreinf $numBarsLayer2 $BarArea 0.0 $steeldistZ 0.0 -$steeldistZ; # 2nd layer reinforcement
	layer straight $IDreinf $numBarsLayer2 $BarArea [expr -$steeldistY/3] $steeldistZ [expr -$steeldistY/3] -$steeldistZ; # 2nd layer reinforcement
	layer straight $IDreinf $numBarsLayer3 $BarArea [expr $steeldistY/3] $steeldistZ [expr $steeldistY/3] -$steeldistZ; # 3rd layer reinforcement
	layer straight $IDreinf $numBarsLayer4 $BarArea  $steeldistY $steeldistZ  $steeldistY -$steeldistZ;				# top layer reinforcement

    };	# end of fibersection definition
puts "Section has been defined"

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
set TransfTag 1; 			# associate a tag to column transformation
set TransfType PDelta;			# options, Linear PDelta Corotational 
geomTransf $TransfType $TransfTag; 	

# element connectivity:
#set eleTag 1;
#set iNode 1;
#set jNode 2;
set numIntgrPts 2;									# number of integration points for force-based element
set memID 1;
set nllength 400.0;

# define elements
for {set x 1} {$x <= $nel } {incr x} {
   set iNode [expr $x]
   set jNode [expr $x + 1]
   element NLDispBeamColumn2d $x $iNode $jNode $numIntgrPts $SecTag $TransfTag;
   puts "$x $iNode $jNode $numIntgrPts $memID $SecTag $TransfTag"
   }


#element nonlinearBeamColumn $eleTag $iNode $jNode $numIntgrPts $SecTag $TransfTag;	# self-explanatory when using variables
puts "Element has been defined"


# Define RECORDERS -------------------------------------------------------------
recorder Node -file DFree.out -time -node $nnode -dof 1 2 3 disp;					# displacements of free nodes
#recorder Node -file $dataDir/DBase.out -time -node 1 -dof 1 2 3 disp;					# displacements of support nodes
recorder Node -file RBase.out -time -node 1 -dof 1 2 3 reaction;				# support reaction
#recorder Element -file $dataDir/FCol.out -time -ele 1 globalForce;					# element forces -- column
#recorder Element -file $dataDir/ForceColSec1.out -time -ele 1 section 1 force;				# Column section forces, axial and moment, node i
#recorder Element -file $dataDir/DefoColSec1.out -time -ele 1 section 1 deformation;			# section deformations, axial and curvature, node i
#recorder Element -file $dataDir/ForceColSec$numIntgrPts.out -time -ele 1 section $numIntgrPts force;	# section forces, axial and moment, node j
#recorder Element -file $dataDir/DefoColSec$numIntgrPts.out -time -ele 1 section 1 deformation;		# section deformations, axial and curvature, node j
#recorder Element -xml $dataDir/PlasticRotation.out -time -ele 1 plasticRotation;			# section deformations, axial and curvature, node j

# # #fiber stress-strain recorders
# set recsec 1;
# recorder Element -file coverstraint.out -time -ele 1 section $recsec fiber $coverY 0.0 $IDconcCover strain
# recorder Element -file coverstresst.out -time -ele 1 section $recsec fiber $coverY 0.0 $IDconcCover stress
# recorder Element -file coverstrainb.out -time -ele 1 section $recsec fiber -$coverY 0.0 $IDconcCover strain
# recorder Element -file coverstressb.out -time -ele 1 section $recsec fiber -$coverY 0.0 $IDconcCover stress

# recorder Element -file corestraint.out -time -ele 1 section $recsec fiber $coreY 0.0 $IDconcCore strain
# recorder Element -file corestresst.out -time -ele 1 section $recsec fiber $coreY 0.0 $IDconcCore stress
# recorder Element -file corestrainb.out -time -ele 1 section $recsec fiber -$coreY 0.0 $IDconcCore strain
# recorder Element -file corestressb.out -time -ele 1 section $recsec fiber -$coreY 0.0 $IDconcCore stress
# recorder Element -file corestrainc.out -time -ele 1 section $recsec fiber 0.0 0.0  $IDconcCore strain
# recorder Element -file corestressc.out -time -ele 1 section $recsec fiber 0.0 0.0  $IDconcCore stress

# recorder Element -file reinfstraint.out -time -ele 1 section $recsec fiber $steeldistY $steeldistZ $IDreinf strain
# recorder Element -file reinfstresst.out -time -ele 1 section $recsec fiber $steeldistY $steeldistZ $IDreinf stress
# recorder Element -file reinfstrainb.out -time -ele 1 section $recsec fiber -$steeldistY $steeldistZ $IDreinf strain
# recorder Element -file reinfstressb.out -time -ele 1 section $recsec fiber -$steeldistY $steeldistZ $IDreinf stress


recorder Element -file profile.out -ele 1 integrationPoints;
 #section deformations, axial and curvature
 set nsec [expr $numIntgrPts*$nel]
 set intincr [expr $numIntgrPts - 1]
for {set x 1} {$x <= $nsec } {incr x $numIntgrPts} {
recorder Element -file DefoColSec$x.out -time -ele [expr ($x+$intincr)/$numIntgrPts] -section 1 deformation;
recorder Element -file DefoColSec[expr $x+1].out -time -ele [expr ($x+$intincr)/$numIntgrPts] -section 2 deformation;
# recorder Element -file DefoColSec[expr $x+2].out -time -ele [expr ($x+$intincr)/$numIntgrPts] -section 3 deformation;
# recorder Element -file DefoColSec[expr $x+3].out -time -ele [expr ($x+$intincr)/$numIntgrPts] -section 4 deformation;
# recorder Element -file DefoColSec[expr $x+4].out -time -ele [expr ($x+$intincr)/$numIntgrPts] -section 5 deformation;
} 		

# define GRAVITY -------------------------------------------------------------
pattern Plain 1 Linear {
   load $nnode 0 $P 0
}

# Gravity-analysis parameters -- load-controlled static analysis
set Tol 1.0e-8;				# convergence tolerance for test
constraints Plain;     			# how it handles boundary conditions
numberer Plain;				# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;			# how to store and solve the system of equations in the analysis
test NormDispIncr $Tol 100; 		# determine if convergence has been achieved at the end of an iteration step
algorithm ModifiedNewton -initial;			# use Newton's solution algorithm: updates tangent stiffness at every iteration
set NstepGravity 20;  			# apply gravity in 10 steps
set DGravity [expr 1./$NstepGravity]; 	# first load increment;
integrator LoadControl $DGravity;	# determine the next time step for an analysis
analysis Static;			# define type of analysis static or transient
analyze $NstepGravity;			# apply gravity

# ------------------------------------------------- maintain constant gravity loads and reset time to zero
loadConst -time 0.0

puts "Model Built"

# STATIC PUSHOVER ANALYSIS --------------------------------------------------------------------------------------------------
# define LATERAL load -------------------------------------------------------------
set F 200e3;
# Lateral load pattern
pattern Plain 2 Linear {;
	load $nnode $F 0.0 0.0;
}

## Analysis parameters
set IDctrlNode $nnode;				# node where displacement is read for displacement control
set IDctrlDOF 1;
set Dincr 0.05;			# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis
set Tol 1.e-5;                 			# Convergence Test: tolerance
set maxNumIter 2000;               		# Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
set printFlag 0;               			# Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
set TestType EnergyIncr;			# Convergence-test type
#set algorithmType ModifiedNewton;

constraints Plain;		
numberer Plain;
system BandGeneral;
test $TestType $Tol $maxNumIter $printFlag;
algorithm ModifiedNewton -initial;        
integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;
analysis Static;

#  ---------------------------------    perform Static Pushover Analysis
set Nsteps [expr 150];		        # number of pushover analysis steps

set currentStep 0;
set ok 0

while {$ok == 0 && $currentStep < $Nsteps} {

	set ok [analyze 1]
	
	if {$ok != 0} {
	
		puts "Trying NewtonWithLineSearch .."
	    algorithm NewtonLineSearch 0.8
	    set ok [analyze 1]
	    if {$ok == 0} {puts "that worked .. back to Broyden"}
	    algorithm ModifiedNewton -initial
	}
        if {$ok != 0} {
	    puts "Trying Newton .."
	    algorithm Newton
	    set ok [analyze 1 ]
	    algorithm ModifiedNewton -initial
	}
	if {$ok != 0} {
	    algorithm Broyden
	    set ok [analyze 1]
	    algorithm ModifiedNewton -initial
	}

	set currentStep [expr $currentStep + 1]
	puts "step $currentStep done"
}

puts "Pushover Done"
wipe;
