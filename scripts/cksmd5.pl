#!/usr/bin/perl
#------------------------------------------------------------------------------#
#                                                                              #
#  cksmd5.pl                                                                   #
#  directory checksum computation perl script                                  #
#                                                                              #
#  (C) SHEIN; Bad Aibling, February 2002                     Steffen Hein      #
#  [ Update: March 23, 2007 ]                             <contact@sfenx.de>   #
#                                                                              #
#------------------------------------------------------------------------------#
#use warnings 'all';
use Env;
use strict;
use cksmd5;

#------------------------------------------------------------------------------#
my ( @files, @args, $directory, $string, $i );

#-----------------------------------------------------------------------------#
# check command line:
#
#if ( ! $ARGV[0] ) {
#    die (
#"Usage: cksmd5.pl DIRECTORY \n"
#    );
#}

#-----------------------------------------------------------------------------#

$directory = "./";
if ( $ARGV[0] ) {
    $directory = $ARGV[0];
}

#------------------------------------------------------------------------------#
# Writing new final package checksum:
#
# check / open  directory :
#
opendir( DIR, $directory )
  || die ( "can't open directory " . $directory . "\n" );
@files = grep -T, readdir(DIR);
close(DIR);

$string = &cksmd5( $directory, "CHECKSUM.MD5" );
print STDOUT ($string);

#-------------------------------------------------------------------------------#
print "Done !\n";
