#   Program:    abYmod gromacs
#   File:       run_gromacs_abymod.pm
#   
#   Version:    V2.0
#   Date:       03.07.20
#   Function:   Main library routines for the abYmod program
#   
#   Copyright:  (c) Prof. Andrew C. R. Martin, UCL, 2013-2020
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#

#use strict;
#use config;
#use util;

my ($inPDB, $tmpDir, $energyOut) = @ARGV;
my $gromacsBinDir= "/usr/local/gromacs/bin";
my $bindir= "/home/lilian/bin";
my $loopflapper="/home/lilian/loop_flapper";


# Gromacs energy minization
my $gromacsDir="$bindir/gromacs";
#my $gromacsParam="gromacsParam";
my $gromacsParamDir="$loopflapper/data/gromacsParams";       # Gromacs parameter files
my $gromacsGenionMdp="$gromacsParamDir/genion.mdp";
my $gromacsEmMdp="$gromacsParamDir/em.mdp";
#my $mmDataDir=


OptimizeModel($inPDB, $tmpDir);
system ("rm -f $tmpDir/run/*");
system ("rm -f $tmpDir/setup/*");

sub RunCommand
{
    my ($exe, $useSystem) = @_;
    my $result = '';

    #print STDERR "$exe\n";
    if(defined($useSystem) && $useSystem)
    {
        $result = system("$exe");
    }
    else
    {
        $result = `$exe`;
    }
    return($result);
}

#>OptimizeModel($inPDB, $tmpDir)
# ------------------------------
# Runs minimization
#
# 28.09.15 Original   By: ACRM
# 02.11.15 Various fixes and enancemanemnts for loopdb
# 09.01.20 Switched out tinker and replaced with gromacs
sub OptimizeModel{
    my($inPDB, $tmpDir) = @_;
    
    # Test for pdb
    #print STDERR "\n\n*** Running Gromacs minimization (this will take a while) ***\n\n";
    if( ! -e $inPDB || ! -s $inPDB ){
        printf STDERR "*** MINIMIZATION FAILED - Input file does not exist ***\n";
        return($inPDB);
    }

    # Setup run
    open(LOG, ">", "$tmpDir/log");
    my $setupDir = "$tmpDir/setup";
    my $runDir = "$tmpDir/run";
    mkdir($setupDir);
    mkdir($runDir);

    # Renumber the PDB file
    my $renumPDB = "$setupDir/renumbered.pdb";
    my $exe0 = "$bindir/pdbrenum $inPDB $renumPDB";
    RunCommand("(cd $setupDir; $exe0)");
    
    # Convert pdb to gromacs
    my $gromacsPDB = "$setupDir/renumbered_gmx.pdb";
    open ($handle, $renumPDB);
    #print $handle
    my $exe1 = "$gromacsBinDir/gmx pdb2gmx -f $renumPDB -water spc -ignh -ff amber03 -o $gromacsPDB 2>&1";
    my $result1 = RunCommand("(cd $setupDir; $exe1)");
    #print STDERR $result1;
    #print LOG $result1;
    if( $result1 !~ /You have successfully generated/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to convert pdb to gromacs ***\n";
        return($inPDB);
    }
   
    # Define solvent box
    my $boxPDB = "$setupDir/box.pdb";
    my $exe2 = "$gromacsBinDir/gmx editconf -f $gromacsPDB -bt cubic -d 0.5 -o $boxPDB 2>&1";
    my $result2 = RunCommand("(cd $setupDir; $exe2)");
    #print STDERR $result2;
    #print LOG $result2;
    if( $result2 !~ /new box vectors/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to define solvent box ***\n";
        return($inPDB);
    }
    
    # Sovate box
    my $solvatedPDB = "$setupDir/solvated.pdb";
    my $exe3 = "$gromacsBinDir/gmx solvate -cp $boxPDB -cs -o $solvatedPDB -p topol.top 2>&1";
    my $result3 = RunCommand("(cd $setupDir; $exe3)");
    #print STDERR $result3;
    #print LOG $result3;
    if( $result3 !~ /Generated solvent/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to solvate box ***\n";
        return($inPDB);
    }
    
    # Prepare to add solvent to topology
    my $exe4 = "$gromacsBinDir/gmx grompp -c solvated.pdb -p topol.top -f $gromacsGenionMdp -o genion.tpr 2>&1";
    my $result4 = RunCommand("(cd $setupDir; $exe4)");
    #print STDERR $result4;
    #print LOG $result4;
    if( $result4 !~ /This run will generate/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to add solvent to topology ***\n";
        return($inPDB);
    }
    
    # Generate neutral solvent system
    my $exe5 = "$gromacsBinDir/gmx genion -s genion.tpr -conc 0.15 -neutral -pname NA -nname CL -o system.gro -p topol.top < $gromacsParamDir/genion_choices.txt 2>&1";
    my $result5 = RunCommand("(cd $setupDir; $exe5)");
    #print STDERR $result5;
    #print LOG $result5;
    if( $result5 !~ /Replacing solvent molecule/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to generate neutral solvent system ***\n";
        return($inPDB);
    }
    
    # Prepare energy minimisation run
    my $exe6 = "$gromacsBinDir/gmx grompp -c  $setupDir/system.gro -p $setupDir/topol.top -f $gromacsEmMdp -o em.tpr 2>&1";
    my $result6 = RunCommand("(cd $runDir; $exe6)");
    #print STDERR $result6;
    #print LOG $result6;
    if( $result6 !~ /This run will generate/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to prepare energy minimisation run ***\n";
        return($inPDB);
    }
        
    # Run minimisation
    print STDERR "\n\n***********\nPlease wait while the energy minimization command below runs. This may take up to 2 minutes depending on your system\n***********\n\n";
    my $minimizedPDB = "$runDir/minimized.pdb";
    my $exe7 = "$gromacsBinDir/gmx mdrun -deffnm em -v -c $minimizedPDB 2>&1";
    my $result7 = RunCommand("(cd $runDir; $exe7)");
    #print STDERR $result7;
    print $result7;

    #open ($handle, '>>', $energyOut);
    #print $handle $result7;
    
    #print LOG $result7;
    if( ! -e $minimizedPDB || $result7 !~ /Energy minimization reached the maximum number of steps before/g ){
        printf STDERR "*** MINIMIZATION FAILED - Failed to run minimization ***\n";
        return($inPDB);
    }
    
}


