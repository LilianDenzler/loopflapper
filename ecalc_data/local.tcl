#*************************************************************************
#
#   Program:    ecalc interface
#   File:       local.tcl
#   
#   Version:    V1.0
#   Date:       29.09.94
#   Function:   Set up local stuff for this installation of ECalc
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1994
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0372) 275775
#   EMail:      INTERNET: amartin@scitec.adsp.sub.org
#                         martin@bsm.bioc.ucl.ac.uk
#               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
#               JANET:    martin@uk.ac.ucl.bioc.bsm
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
#*************************************************************************
#
#   Revision History:
#   =================
#   V0.1  29.09.94 Original    By: ACRM
#   V1.0  29.09.94 First release version
#
#*************************************************************************
# Set the ecalc path for running the main program
# -----------------------------------------------
set ecalc "~amartin/bin/ecalc"


##########################################################################
# Add the ecalc directory to the execute path
# -------------------------------------------
set auto_path "~amartin/ecalc $auto_path"


##########################################################################
# Define colours for the Run and Quit buttons. Change these for B&W
# screens!
set ActiveColour  Red
set PassiveColour Cyan
set AdvColour     GreenYellow


##########################################################################
# Set maximum number of disulphides, zones and ignores
set maxss      12
set maxzones    8
set maxignores  8
