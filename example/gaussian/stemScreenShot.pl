#!/usr/bin/perl
$filename = $ARGV[0];
$xrotation = $ARGV[1];
$zrotation = $ARGV[2];
$scale = $ARGV[3];
$topView = $ARGV[4];

$output = $filename;
$filename = $filename . ".zcon";
$output = $output . "-";
$output = $output . $xrotation;
$output = $output . "-". $zrotation . ".gif";
open my $fh, '>', 'gnuplot.in' or die "Cannot open file: $!";
if(defined($topView))
{
print $fh "set pm3d map\n";
print $fh "set terminal gif size 800, 800\n";
print $fh "set size square\n";
print $fh "set output \"$output\"\n";
print $fh "set palette defined (0 \"black\", 10 \"white\")\n";
print $fh "splot \"$filename\"\n";
}
else
{
print $fh "set pm3d\n";
print $fh "set hidden3d\n";
print $fh "set term gif\n";
print $fh "set size square\n";
print $fh "set zrange[10:5000]\n";
print $fh "set output \"$output\"\n";
print $fh "set view $xrotation, $yrotation, $scale \n";
print $fh "splot \"$filename\" with lines\n";
}
close $fh;
system "gnuplot gnuplot.in";
