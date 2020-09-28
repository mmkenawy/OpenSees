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

set nel 6
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
 	

section Fiber $SecTag {;		# Define the fiber section
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

# define elements
for {set x 1} {$x <= $nel } {incr x} {
   set iNode [expr $x]
   set jNode [expr $x + 1]
   element NLDispBeamColumn2d $x $iNode $jNode $numIntgrPts $SecTag $TransfTag 1 -nllength 400.0;
   puts "$x $iNode $jNode $numIntgrPts $SecTag $TransfTag"
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
#recorder Element -file DefoColSec[expr $x+2].out -time -ele [expr ($x+$intincr)/$numIntgrPts] -section 3 deformation;
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
#set Dincr 0.1;			# displacement increment for pushover. you want this to be very small, but not too small to slow down the analysis
set Tol 1.e-4;                			# Convergence Test: tolerance
set maxNumIter 2000;               		# Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
set printFlag 0;               			# Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
set TestType EnergyIncr;			# Convergence-test type
#set algorithmType ModifiedNewton;

constraints Plain;		
numberer Plain;
system BandGeneral;
test $TestType $Tol $maxNumIter $printFlag;
#algorithm ModifiedNewton -initial;        
#integrator DisplacementControl  $IDctrlNode $IDctrlDOF $Dincr;

set incrnum 0;

foreach Dincr {-0.00093127  0.00093127   0.0055877   0.0018624   0.0027939    0.005588    0.003725    0.005587    0.007452    0.003723    0.005588     0.00745    0.008382   -0.001862   -0.003726   -0.005588   -0.002794   -0.002794   -0.002794   -0.002794   -0.004656   -0.003723   -0.004658   -0.003724   -0.004657   -0.003725  -0.0018625  -0.0065189  -0.0027939  -0.0055876  -0.0018626  -0.0027935   -0.001863   -0.003725   -0.005587    -0.00652   -0.005588   -0.003723   -0.007453   -0.007449   -0.003724   -0.002794   -0.005588   -0.003726    0.001864    0.004656    0.004656    0.004656    0.004658    0.002791    0.008382    0.005588    0.005588     0.00745   0.0093133    0.006519   0.0027937   0.0074501   0.0083819    0.006518    0.010246     0.00745    0.012105    0.008382    0.010244    0.012106    0.009314    0.013033     0.01118     0.01583     0.01397     0.01583     0.00652    -0.00093    -0.00186     -0.0028    -0.00745    -0.00838    -0.01769    -0.01956   -0.015831    -0.01397   -0.012106   -0.012106   -0.011176   -0.008382   -0.019556  -0.0093129   -0.012106  -0.0009314  -0.0093123   -0.004657   -0.007451   -0.008379   -0.009314   -0.008382    -0.00838   -0.009314   -0.009312   -0.012108   -0.017693    -0.01211    -0.00838    -0.01304     -0.0121    -0.00932    -0.00558    -0.00187      0.0028     0.00559     0.01396     0.00652      0.0149     0.01211     0.01211    0.013035    0.009312    0.008382    0.012108     0.00745    0.011173    0.008382    0.009314    0.006518    0.005588    0.010244    0.015832    0.016764    0.013038    0.013037    0.019556    0.010244    0.010244    0.013969     0.01397     0.01956     0.01583     0.00745    -0.00186    -0.00652    -0.00838    -0.00932    -0.02887   -0.016759    -0.01397    -0.01397   -0.017694   -0.027007   -0.017694  -0.0018626  -0.0046563  -0.0065189   -0.011176   -0.010245   -0.009311   -0.010244   -0.008382   -0.009312    -0.02887   -0.009309    -0.01025    -0.01304    -0.01117    -0.01118    -0.01397    -0.00372    -0.00093     0.01304     0.01397     0.01396      0.0177    0.017693    0.023282    0.015832    0.011173    0.018627    0.016763    0.012106    0.011176    0.027008    0.009311      0.0149     0.01397    0.010244    0.017695     0.01583      0.0149     0.01304     0.01304     0.01676     0.01956     0.01862     0.01956     0.03631     0.03167      0.0298     0.00279     0.00094    -0.00465    -0.01676    -0.00747    -0.00929    -0.01118    -0.01026    -0.01489    -0.01397    -0.01304    -0.01863      -0.027    -0.00745    -0.01397    -0.01304    -0.01397    -0.01769    -0.01863   -0.013037   -0.015832   -0.016764   -0.016762   -0.014899    -0.02794   -0.013969   -0.010244   -0.020489   -0.015831   -0.013968   -0.017696   -0.013039    -0.01676    -0.01304    -0.01304    -0.01303    -0.01584    -0.01769    -0.01956    -0.02421    -0.08195    -0.00932    -0.01117    -0.00559    -0.00371     0.00279      0.0093      0.0177     0.02142      0.0242     0.02702     0.02328     0.02514     0.02887     0.04377    0.028873    0.030732    0.027008    0.029799     0.02887    0.025143    0.026076    0.026079     0.02235     0.03073     0.02421     0.02701     0.02514       0.027     0.03726       0.027     0.01585     0.00186    -0.00744    -0.00747    -0.00744    -0.01024    -0.01676    -0.01585    -0.06798    -0.04936    -0.04842   -0.019559   -0.018626   -0.023282   -0.031663   -0.016763   -0.016763   -0.039113   -0.024214   -0.031667      -0.027    -0.02608    -0.04098    -0.02607    -0.02142    -0.02979    -0.03353      -0.027     0.00371     0.00465     0.01117     0.01397     0.01585     0.02235     0.01768     0.01956     0.01863     0.02142     0.02142     0.02328     0.04936    0.027939    0.026076    0.021417    0.027008    0.025144    0.024213     0.03632    0.039113     0.02608     0.02607     0.02701      0.0298     0.03259     0.02794      0.0298     0.03073      0.0298     0.02235     0.03261     0.02606      0.0298     0.00467     0.00091     0.00188    -0.00467    -0.00371    -0.00653    -0.00652    -0.02048    -0.02888    -0.01303    -0.03726    -0.01488    -0.02144    -0.02047    -0.02141    -0.02236    -0.02329    -0.01769    -0.02608    -0.02328     -0.0298    -0.02421    -0.02794   -0.027007   -0.040978   -0.022349   -0.015832   -0.062395   -0.033526   -0.021423    -0.02328    -0.02049    -0.02514    -0.02049    -0.01955    -0.02981    -0.02421    -0.01582    -0.01585      -0.027    -0.03074    -0.03723    -0.03074    -0.02329    -0.02886    -0.00746    -0.01024    -0.00653      0.0028     0.00838     0.01491     0.01676     0.01768     0.03632     0.02982     0.02886      0.0242      0.0233     0.03167     0.02514     0.02886     0.02981     0.05308    0.042837    0.060534    0.021419     0.02235    0.031663    0.029802    0.049355     0.03912       0.027      0.0298     0.04471     0.04191     0.02047     0.03168     0.02979     0.02421     0.03911     0.03074     0.02514    -0.00094    -0.00279    -0.01676    -0.00836    -0.01491    -0.01491    -0.01303    -0.01582    -0.01585    -0.01303    -0.01677     -0.0177    -0.01862    -0.02606    -0.02329    -0.01397    -0.01491    -0.01768    -0.02608    -0.01676    -0.01863    -0.02608    -0.02235    -0.02328   -0.027937   -0.019556    -0.02049   -0.012105   -0.056807   -0.023282    -0.02142   -0.024213    -0.05215    -0.04098    -0.02235    -0.02794    -0.02793     -0.0354    -0.03165    -0.04565    -0.04376    -0.02235    -0.02327     -0.0177    -0.00838    -0.00745    -0.00094     0.00465     0.01118     0.01676     0.01583     0.04937     0.01862       0.027     0.01865      0.0307     0.02424     0.03073     0.03724     0.03073     0.04657     0.04004    0.021421    0.064259    0.067051    0.031664    0.027935      0.0326     0.03632     0.03725     0.03725     0.04192     0.03352     0.02421     0.03165     0.02982     0.04562     0.03446     0.04471     0.03723     0.03912     0.03818     0.02794     0.00279     0.00094    -0.00094    -0.01488    -0.02329    -0.03259    -0.02609    -0.01676    -0.02236    -0.01861    -0.02888    -0.02792    -0.02235    -0.02609    -0.02794      -0.054    -0.07452    -0.04376    -0.08475    -0.03166   -0.049357   -0.027939   -0.025144   -0.050289   -0.063331    -0.05401    -0.02887    -0.03352    -0.09035    -0.06982    -0.06056    -0.02047    -0.03538    -0.07638    -0.04003    -0.03912    -0.01026    -0.00838    -0.00559    -0.00091     0.00185     0.01303     0.01397     0.01397     0.01956     0.02608     0.02886     0.04379      0.0717     0.04562     0.04006     0.03817     0.06053     0.03353     0.05495     0.03446    0.059599    0.081022     0.13876     0.08847     0.04099     0.04003     0.04285     0.03911     0.04003     0.03912     0.04376     0.03912     0.05867     0.02329     0.03259     0.02142    -0.00651    -0.00652    -0.01491     -0.0298    -0.02235    -0.03818    -0.04097    -0.01864    -0.02421    -0.02047    -0.02144    -0.01676    -0.01862    -0.01956     -0.0177    -0.01956    -0.02327     -0.0177    -0.02047    -0.03168    -0.02328    -0.02421    -0.01676    -0.02701   -0.083815   -0.079158    -0.05215   -0.048427    -0.04191    -0.04191    -0.04935    -0.03538    -0.03726    -0.06985     -0.0503    -0.07635    -0.04376      -0.068    -0.02979    -0.01677    -0.02049    -0.00186      0.0028     0.01117     0.02329     0.01677     0.03165     0.03446     0.03259     0.03912     0.03353     0.04749     0.06147     0.02886     0.03817     0.04844     0.05494    0.090337     0.22164     0.11175     0.06893     0.05679     0.05682     0.04471      0.0717     0.04283     0.05308     0.05494     0.04285     0.03632     0.03818     0.02141     0.00932     0.00559     0.00094     0.00092    -0.00186    -0.01582    -0.02888    -0.01862    -0.01491    -0.01956    -0.03073    -0.04562     -0.0475    -0.03538    -0.05123    -0.04656    -0.03632    -0.04006    -0.16017    -0.09872   -0.063323   -0.013038   -0.012106   -0.010244   -0.003726   -0.008379   -0.005588   -0.010244   -0.012108   -0.015831   -0.011175   -0.016763      -0.149    -0.05122    -0.01397    -0.03446    -0.01491    -0.01303    -0.01861     -0.0205    -0.01397    -0.02047    -0.07453    -0.03071    -0.02143    -0.01583    -0.03353    -0.03817    -0.30361    -0.09312
} {
	set incrnum [expr $incrnum + 1];
	puts "incr $incrnum";

  algorithm ModifiedNewton -initial
  #algorithm NewtonLineSearch
  integrator DisplacementControl $IDctrlNode 1 $Dincr
  analysis Static;
  #analyze 100
  
set Nsteps 100;		        # number of steps per incr

set currentStep 0;
set ok 0

while {$ok == 0 && $currentStep < $Nsteps} {

	set ok [analyze 1]
	
	if {$ok != 0} {
	
		puts "Trying NewtonWithLineSearch .."
	    algorithm NewtonLineSearch 0.8
		#algorithm ModifiedNewton -initial
		test $TestType $Tol $maxNumIter $printFlag;
	    set ok [analyze 1]
	    if {$ok == 0} {puts "that worked .. back to Broyden"}
	    algorithm ModifiedNewton -initial
		#algorithm NewtonLineSearch
	}
        if {$ok != 0} {
	    puts "Trying Newton .."
	    algorithm Newton
		test $TestType $Tol $maxNumIter $printFlag;
	    set ok [analyze 1 ]
	    algorithm ModifiedNewton -initial
		#algorithm NewtonLineSearch
	}
	if {$ok != 0} {
	    algorithm Broyden
		#algorithm ModifiedNewton -initial
		test $TestType $Tol $maxNumIter $printFlag;
	    set ok [analyze 1]
	    algorithm ModifiedNewton -initial
		#algorithm NewtonLineSearch
	}

	set currentStep [expr $currentStep + 1]
#	puts "step $currentStep done"
}
}
puts "Cyclic test Done"
wipe;