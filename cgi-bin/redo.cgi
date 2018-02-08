#!/bin/perl

# -----------------------------------------------------------------------
#  File Name: redo.cgi  v.8   
#
#  Description: This Perl CGI removes the input filename 
#               from the Phred Output Data Directory (PODD).  
#               
#  Created: 04/07/2000
#  Modified: 05/17/2000
#
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib";
&Parse_Form;

$file = $formdata{'file'};
$dir = "$FASTA_DIR";

$fasta = "$dir/$file";
$qual  = "$dir/$file.qual";
$hist  = "$dir/$file.hist";

system "rm -f $fasta $qual $hist";

sleep 5;
print "Location: fastadisplay.cgi\n\n";

# -----------------------------------------------------------------------
# For more information see:
# Lazo, G.R., Tong, J., Miller, R., Hsia, C., Rausch, C., Kang, Y., 
# and Anderson, O.D. 2001. Software scripts for quality checking of 
# high-throughput nucleic acid sequencers. BioTechniques Vol 30, No. 6.
#
# This software is a "United States Government Work" under the terms 
# of the United States Copyright Act. It was written as part of official 
# duties performed by a United States Government employee(s) and thus 
# cannot be copyrighted. This software is freely available to the public 
# for use. The United States Department of Agriculture and the U.S. 
# Government have not placed any restriction on its use or reproduction. 
# -----------------------------------------------------------------------
