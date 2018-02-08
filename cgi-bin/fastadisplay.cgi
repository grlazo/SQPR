#!/bin/perl

# -----------------------------------------------------------------------
#  File Name: fastadisplay.cgi v.8                              
#                                                       
#  Description: This Perl cgi-bin script lists links to  
#               the Phred Output Data Directory (PODD). 
#               Selecting an entry initiates the phddisplay.cgi 
#               script to display statistics for the selected run. 
#                                                       
#  Created: 02/22/2000
#  Modified: 07/12/2000
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib";
&Parse_Form;

print "Content-type: text/html\n\n";

if (exists $formdata{'dir'}) {
   $file = $formdata{'dir'};
   $frmscript = $formdata{'script'};
}

$dir = $file . '.fasta';
$origin = $ENV{'HTTP_REFERER'};

print <<"HEADER HTML";
<HTML>
<HEAD>
<TITLE>SQPR: PHRED-Analyzed Sequencing Runs</TITLE>
</HEAD>

<BODY BGCOLOR="FFFFCC">
<H1>SQPR:</H1>
<H2>PHRED-Analyzed Sequencing Runs</H2>
HEADER HTML

print <<"HEAD2 HTML";
<HR>
<TABLE BORDER=0 WIDTH=600> 
<TR><TD ALIGN=RIGHT>
<A HREF=\"$CGI_BIN/dirdisplay.cgi\" style=\"text-decoration:none\">Click here if you can't find your file</A>
</TD></TR></TABLE>
<HR>
HEAD2 HTML

chdir ("$FASTA_DIR");

opendir (DIR, ".") || &ErrorMessage;
@files = readdir (DIR);
closedir (DIR);

@files = sort{$b cmp $a}(@files);

if (@files) {
   print "<P><H3>Select a file to display:</H3><BR>";
   if ($frmscript eq "$CGI_BIN/runphred.cgi") {
      foreach $filename (@files) {
	 if ($filename =~ /fasta$/) {
	   if ($filename eq $dir) {
		print "<LI><A HREF=\"$CGI_BIN/phddisplay.cgi?file=$filename\" style=\"text-decoration:none\">$filename</A> <-- This file was just processed.<BR>\n" unless ($filename =~ /^\.+$/);
	   } else {
	        print "<LI><A HREF=\"$CGI_BIN/phddisplay.cgi?file=$filename\" style=\"text-decoration:none\">$filename</A><BR>\n" unless ($filename =~ /^\.+$/);
	   }
         }
      }
   } else {
      foreach $filename (@files) {
         if ($filename =~ /fasta$/) {
           print "<LI><A HREF=\"$CGI_BIN/phddisplay.cgi?file=$filename\" style=\"text-decoration:none\">$filename</A><BR>\n" unless ($filename =~ /^\.+$/);
         }
      }
   }
}

print <<"FOOTER HTML";
</BODY>
</HTML>
FOOTER HTML

sub ErrorMessage {
	print "Content-type: text/html\n\n";
	print "The server can't access the file. Aborting script. \n";
	exit;
}

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
