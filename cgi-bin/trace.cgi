#!/bin/perl

# -----------------------------------------------------------------------
# Filename: trace.cgi v.8
#
# Description: This Perl cgi-bin script takes in the 
#              directory name and the trace file name.  
#              With the selected trace file data, a 
#              viewable trace file is created.
#
# Created: 04/05/2000
# Modified: 04/11/2000
# -----------------------------------------------------------------------

require "GlobalVariables.lib";
require "subparseform.lib"; 
&Parse_Form;

print "Content-type: text/html\n\n";

$directory = $formdata{'filename'};
$well = $formdata{'well'};

$dir = "$TRACE_DIR/$directory";

$welloc = &deWELL($well);					# Decode well location
&findtrace($directory,$welloc);					# find trace file with welloc

$fname = "$WEB_TRACE_LINK_DIR/$tracefile";

system "rm -f $WEB_TRACE_DIR/*";				# Clear all the links in the www dir
system "ln -s $dir/$tracefile $WEB_TRACE_DIR/$tracefile";	# create a link to a www dir

print <<"HTML code";
<HTML>
<HEAD>
<TITLE> $tracefile </TITLE>
</HEAD>
<BODY BGCOLOR="FFFFCC"> 
<CENTER> $tracefile
<P>

<applet codebase="$JAVA_APPLET_CLASS_URL" code="ChromatogramApplet.class" width=450 height=200>
<param name="file" value="$fname">
</applet>
<br><font size=-1>
Chromatogram Applet, Release 1, 6/30/96 by Eugen Buehler
</font>

</CENTER>
</BODY>
</HTML>
HTML code

#---------Subroutines---------------------------

sub deWELL {
        $w = $_[0];

        if (($w >= 1) && ($w <= 12))    {
           if ($w <= 9) {$value = '0'.$well; $welloc = 'A'.$value;}
           else {$value = $well; $welloc = 'A'.$value;}
        } elsif (($w >= 13) && ($w <= 24)) {
           $v = $well - 12;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'B'.$value;}       
           else {$value = $v; $welloc = 'B'.$value;}        
        } elsif (($w >= 25) && ($w <= 36)) {
           $v = $well - 24;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'C'.$value;}       
           else {$value = $v; $welloc = 'C'.$value;}        
        } elsif (($w >= 37) && ($w <= 48)) {
           $v = $well - 36;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'D'.$value;}       
           else {$value = $v; $welloc = 'D'.$value;}        
        } elsif (($w >= 49) && ($w <= 60)) {
           $v = $well - 48;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'E'.$value;}       
           else {$value = $v; $welloc = 'E'.$value;}        
        } elsif (($w >= 61) && ($w <= 72)) {
           $v = $well - 60;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'F'.$value;}       
           else {$value = $v; $welloc = 'F'.$value;}        
        } elsif (($w >= 73) && ($w <= 84)) {
           $v = $well - 72;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'G'.$value;}       
           else {$value = $v; $welloc = 'G'.$value;}        
        } elsif (($w >= 85) && ($w <= 96)) {
           $v = $well - 84;             
           if ($v <= 9) {$value = '0'.$v; $welloc = 'H'.$value;}       
           else {$value = $v; $welloc = 'H'.$value;}  
	  }      
}

sub findtrace {
	$d = $_[0];
	$w = $_[1];

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
                if ($w eq $trwell) {$tracefile = $trace;}        # get the trace file name
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

