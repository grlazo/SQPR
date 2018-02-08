#!/bin/perl

# -----------------------------------------------------------------------
#  Filename: trace2.cgi v.8
#
#  Description: This Perl CGI script takes in the directory 
#		name and the trace file name.  With that 
#		trace file data, it will create a viewable 
#		trace file for the 384-format (8x12) plate.
#
#  Created: 08/01/2000
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib"; 
&Parse_Form;

print "Content-type: text/html\n\n";

$directory = $formdata{'filename'};
$well = $formdata{'well'};
$startn = $formdata{'startn'};
$startl = $formdata{'startl'};
#print "well = $well<BR>";
#print "startn = $startn<BR>";
#print "startl = $startl<BR>";

$dir = "$TRACE_DIR/$directory";

@setl1 = ("A", "C", "E", "G", "I", "K", "M", "O");
@setl2 = ("B", "D", "F", "H", "J", "L", "N", "P");
@setn1 = ("01","03","05","07","09","11","13","15","17","19","21","23");
@setn2 = ("02","04","06","08","10","12","14","16","18","20","22","24");

$welloc = &deWELL($well,$startn,$startl);		# Decode well location
#print "directory = $directory<BR>";
#print "welloc = $welloc<BR>";

$tracefile = &findtrace($directory,$welloc);		# find trace file with welloc
#print "trace = $tracefile<BR>";

$fname = "$WEB_TRACE_LINK_DIR/$tracefile";

system "rm -f $WEB_TRACE_DIR/*";		# Clear all the link in the www dir
system "ln -s $dir/$tracefile $WEB_TRACE_DIR/$tracefile";	#create a link to a www dir

print <<"HTML code";
<HTML>
<HEAD>
<TITLE> $tracefile </TITLE>
</HEAD>
<BODY BGCOLOR="FFFFCC"> 
<CENTER> $tracefile
<P>

<applet codebase="$JAVA_APPLET_CLASS_URL" code="ChromatogramApplet.class" width=400 height=200>
<param name="file" value="$fname">
</applet>

</CENTER>
</BODY>
</HTML>
HTML code


#---------Subroutines---------------------------
sub deWELL {
        $w = $_[0];		# well index
	$n = $_[1];		# starting well number for 384-format
	$l = $_[2];		# starting well letter
	

      if ($l eq 'A') {@setletter = @setl1;}
      elsif ($l eq 'B') {@setletter = @setl2;}

      if ($n == 1) {@setnum = @setn1;}
      elsif ($n == 2) {@setnum = @setn2;}

	$r = $w%12;		# remainer
	$d = $w/12;		# dividend

	if ($r == 0) {
		$indexnum = 11;
		$num = $setnum[0];
		$let = $setletter[$indexnum];
		return $welloc = $let.$num;
	} else {
		$indexnum = $d;
		$num = $setnum[$r-1];
		$let = $setletter[$indexnum];
		return $welloc = $let.$num;
	} 
}

sub findtrace {
	$d = $_[0];		# directory name
	$w = $_[1];		# well location

        chdir ("$TRACE_DIR/$d");       
        
        opendir(TR, ".") || &ErrorMessage;
        @trfiles = readdir(TR);
        closedir(TR);
        shift(@trfiles); shift(@trfiles);
        if (@trfiles) {
           foreach $eatr (@trfiles) {
		$trace = $eatr;
                $eatr =~ tr/a-z/A-Z/;
                @nametr = split(/\_/, $eatr);
                $trwell = $nametr[1];
                if ($w eq $trwell) { return $tracefile = $trace;}        # get the trace file name
           }
        }   
}

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

