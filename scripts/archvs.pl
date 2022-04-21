#!/usr/bin/perl
#------------------------------------------------------------------------------#
#                                                                              #
#  new-arch.pl                                                                 #
#  elfe archive creation perl script                                           #
#                                                                              #
#  (C) SHEIN; Bad Aibling, November 2007                     Steffen Hein      #
#  [ Update: March 30, 2022 ]                             <contact@sfenx.de>   #
#                                                                              #
#------------------------------------------------------------------------------#
#use warnings 'all';
BEGIN { push @INC, './'}
use Env;
use strict;
use cksmd5;
$ENV{OSTYPE}="linux";

#------------------------------------------------------------------------------#
my ( @dirtree,  @path,    @files,    @args );
my ( @tarlist,  @srclist, @ziplist,  @tgzlist, @modlist, @text );
my ( $portname, $release, $distname, $directory, $infile, $outfile, $cksmd5 );
my ( $srcarch, $tararch, $tgzarch, $modarch, $gzarch, $bz2arch, $i, $j, $md5 );

#------------------------------------------------------------------------------#
$portname = "danse";
$release  = "1.0r3";
$distname = $portname . "-" . $release;
$srcarch  = $distname . ".src.tar";
$tararch  = $distname . ".tar";
$tgzarch  = $distname . ".tgz";
$gzarch   = $srcarch . ".gz";
$bz2arch  = $srcarch . ".bz2";
$outfile  = "CHECKSUM.MD5";

# for OSTYPE=freebsd, e.g.
# $md5=0/1/2: [ don't use / use only / use also ]
# FreeBSD's implementation of checksum function md5 [ for OSTYPE=freebsd ]
# GNU/Linux implementation of checksum function md5sum [ for OSTYPE=linux ]
$md5 = 2;

#------------------------------------------------------------------------------#
# the directory tree of $distname 
#
$i = 0;
$dirtree[ $i++ ] = $distname;
$dirtree[ $i++ ] = $distname . "/bin";
$dirtree[ $i++ ] = $distname . "/doc";
$dirtree[ $i++ ] = $distname . "/solver";
$dirtree[ $i++ ] = $distname . "/former";
$dirtree[ $i++ ] = $distname . "/poster";
$dirtree[ $i++ ] = $distname . "/math";
$dirtree[ $i++ ] = $distname . "/models";
$dirtree[ $i++ ] = $distname . "/objects";
$dirtree[ $i++ ] = $distname . "/samples";
$dirtree[ $i++ ] = $distname . "/scripts";
$dirtree[ $i++ ] = $distname . "/tools";
$dirtree[ $i++ ] = $distname . "/work";

#------------------------------------------------------------------------------#
# the list of files stored in $tgzarch
#
$i = 0;
$tgzlist[ $i++ ] = "README";
$tgzlist[ $i++ ] = "LICENSE";
$tgzlist[ $i++ ] = "INSTALL";
	#$tgzlist[ $i++ ] = "Makefile";
$tgzlist[ $i++ ] = $distname . "/SETUP";
$tgzlist[ $i++ ] = $distname . "/INSTALL";
$tgzlist[ $i++ ] = $distname . "/LICENSE";
$tgzlist[ $i++ ] = $distname . "/README";
$tgzlist[ $i++ ] = $distname . "/CONFIG.H";
$tgzlist[ $i++ ] = $distname . "/bin";
$tgzlist[ $i++ ] = $distname . "/doc";
$tgzlist[ $i++ ] = $distname . "/solver";
$tgzlist[ $i++ ] = $distname . "/former";
$tgzlist[ $i++ ] = $distname . "/poster";
$tgzlist[ $i++ ] = $distname . "/math";
$tgzlist[ $i++ ] = $distname . "/scripts";
$tgzlist[ $i++ ] = $distname . "/tools";
$tgzlist[ $i++ ] = $distname . "/Makefile";
$tgzlist[ $i++ ] = $distname . "/makefile.unx";
$tgzlist[ $i++ ] = $distname . "/mk.unx";
$tgzlist[ $i++ ] = $distname . "/model.c";
$tgzlist[ $i++ ] = $distname . "/models";
$tgzlist[ $i++ ] = $distname . "/objects/.directory";
$tgzlist[ $i++ ] = $distname . "/samples";
$tgzlist[ $i++ ] = $distname . "/work";

#------------------------------------------------------------------------------#
# the DSC model archives $modlist;
#
#$i = 0;
#$modlist[ $i++ ] = $distname . "/models/mod1";
#$modlist[ $i++ ] = $distname . "/models/mod2";
#$modlist[ $i++ ] = $distname . "/models/mod3";
#$modlist[ $i++ ] = $distname . "/models/mod4";
#$modlist[ $i++ ] = $distname . "/models/mod5";
#$modlist[ $i++ ] = $distname . "/models/mod6";

#------------------------------------------------------------------------------#
# The Z-compressed archive of the following [ essential ] sources; $tararch:
#
$i = 0;
$ziplist[ $i++ ] = "README";
$ziplist[ $i++ ] = "LICENSE";
$ziplist[ $i++ ] = "INSTALL";
	#$ziplist[ $i++ ] = "Makefile";
$ziplist[ $i++ ] = $distname . "/SETUP";
$ziplist[ $i++ ] = $distname . "/INSTALL";
$ziplist[ $i++ ] = $distname . "/LICENSE";
$ziplist[ $i++ ] = $distname . "/README";
$ziplist[ $i++ ] = $distname . "/CONFIG.H";
$ziplist[ $i++ ] = $distname . "/bin";
$ziplist[ $i++ ] = $distname . "/doc";
$ziplist[ $i++ ] = $distname . "/solver";
$ziplist[ $i++ ] = $distname . "/former";
$ziplist[ $i++ ] = $distname . "/poster";
$ziplist[ $i++ ] = $distname . "/math";
$ziplist[ $i++ ] = $distname . "/scripts";
$ziplist[ $i++ ] = $distname . "/tools";
$ziplist[ $i++ ] = $distname . "/Makefile";
$ziplist[ $i++ ] = $distname . "/makefile.unx";
#$ziplist[ $i++ ] = $distname . "/makefile.lnx";
$ziplist[ $i++ ] = $distname . "/mk.unx";
	#$ziplist[ $i++ ] = $distname . "/mk.lnx";
$ziplist[ $i++ ] = $distname . "/model.c";
$ziplist[ $i++ ] = $distname . "/models";
	#$ziplist[ $i++ ] = $distname . "/models/model.c";
$ziplist[ $i++ ] = $distname . "/objects/.directory";
$ziplist[ $i++ ] = $distname . "/samples";

#------------------------------------------------------------------------------#
# Create a new tar archive for FreeBSD port building system:
#
$i = 0;
$srclist[ $i++ ] = "README";
$srclist[ $i++ ] = "LICENSE";
$srclist[ $i++ ] = "INSTALL";
	#$srclist[ $i++ ] = "Makefile";
$srclist[ $i++ ] = $distname;
#------------------------------------------------------------------------------#
# Writing special directory checksums:
#
print "\nwriting directory checksums:\n";
if ( $md5 ) {

    $outfile = "CHECKSUM.MD5";

    if ( ($md5) && ( $ENV{OSTYPE} =~ m/freebsd/ ) ) {
       @args = ("./cksmd5.sh $outfile");
    } else {
       @args = ("./md5sum.sh $outfile");
    }
    system(@args);
}    # end if

# switch back to package directory, then write directory [ tree ] checksums:
#
chdir "../../";
if ( $md5 != 1 )
{
    $i = 0;
    while ( $dirtree[$i] ) {

        # check / open  directory :
        $directory = $dirtree[$i];

        if ( opendir( DIR, $directory ) ) {
            @files = grep -T, readdir(DIR);
            close(DIR);

            @path = split ( /\//, $directory );

            $j = 0;
            until ( !$path[$j] ) {
                $outfile = $path[$j] . ".md5";
                $j++;
            }

            #         unlink $directory ."/". $outfile;
            #         unlink $directory ."/*.MD5";
            #         unlink $directory ."/*.md5";

            $cksmd5 = &cksmd5( $directory, $outfile );

            open( FILE, ">" . $directory . "/" . $outfile );
            print FILE ($cksmd5);
            close(FILE);
        }
        else {
            print "\ncan't open directory " . $directory;
        }
        $i++;
    }    # end while
}    # end if $ENV ...
#------------------------------------------------------------------------------#
# Creating achive for selected files: 
#
# switch to script directory then create archive:
#
print "\ncreating archives:\n";
print "tar -czf spcf-". $release. ".tgz\n";
chdir $distname . "/scripts";
@args = ("./new-spcf.sh");
system(@args);

#------------------------------------------------------------------------------#
# Create a new tgz archive [ program package ], $tgzarch:
#
# switch back to package directory, then create archives:
#
chdir "../../";

#print "\ncreating tgz-archive:";
print "\ntar -czf " . $tgzarch;
@args = ("tar -czf $tgzarch @tgzlist");
system(@args);

#print "\narchive ".$tgzarch." ready !";
#------------------------------------------------------------------------------#
# Create a new Z-compressed archive of essential sources [ $tararch.".Z" ]:
#
# switch into package directory, then create archive:
print "\ncreating archive ".$tararch.".Z:";

print "\ntar -cf " . $tararch;
@args = ("tar -cf $tararch @ziplist");
system(@args);
print "\ncompressing to ~.Z";
@args = ("compress -f $tararch");
system(@args);

print "\narchive ".$tararch.".Z ready !";
#------------------------------------------------------------------------------#
# Create a new tar archive for FreeBSD port building system:
#
#print "\ncreating tar archives ".$srcarch.":";
#print "\ncreating src/release
#tar archives:";

print "\ntar -cf " . $srcarch;
@args = ("tar -cf $srcarch @srclist");
system(@args);
print "\ngzipping into ~.gz file";
@args = ("gzip -cf $srcarch > $gzarch");
system(@args);

#print "\narchive ".$gzarch." ready !";
print "\nbzipping into ~.bz2 file";
@args = ("bzip2 -f $srcarch");
system(@args);

#print "\narchive ".$srcarch.".bz2 ready !";
#------------------------------------------------------------------------------#
# model archives:
#
#print "\n\ncreating DSC model archives:";
#$i = 0;
#$j = 1;
#
#while ( $modlist[$i] ) {
#    $modarch = "mod$j-" . $release . "\.tgz";
#    print "\n". $modlist[$i] . " -> ";
#    print $modarch;
#    @args = ("tar -czf $modarch $modlist[$i]");
#    system(@args);
#
#    print "\n archive ".$modarch." ready !";
#    $i++;
#    $j++;
#}
#
#------------------------------------------------------------------------------#
# Writing new final package checksum:
#
$directory = "./";
print "\n\nwriting final package checksum\n";

if ( $md5 ) {

    $outfile = "CHECKSUM.MD5";
    @args = ("rm -f $outfile");
    system(@args);

    if ( ($md5) && ( $ENV{OSTYPE} =~m/freebsd/ ) ) {
       @args = ("md5 * > $outfile");
    } else {
       @args = ("md5sum * > $outfile");
    }
    system(@args);
}

# check / open  directory :
#
# outcommented 05.09.2002:
# if ( ( ( $md5 != 1 ) && ( $ENV{OSTYPE} =~ m/freebsd/ ) )
#   || ( $ENV{OSTYPE} =~ m/linux/ ) )
if ( $md5 != 1 ) # since 05.09.2002
{
    $outfile = $directory . ".md5";

    #  unlink $directory ."/". $outfile;
    #  unlink $directory ."/*.MD5";
    #  unlink $directory ."/*.md5";

    if ( opendir( DIR, $directory ) ) {
        @files = grep -T, readdir(DIR);
        close(DIR);

        $cksmd5 = &cksmd5( $directory, $outfile );

        open( FILE, ">" . $directory . "/" . $outfile );
        print FILE ($cksmd5);
        close(FILE);
    }
    else {
        print "\ncan't open directory " . $directory;
    }    # end if ( opendir )
}    # end if $ENV...

#------------------------------------------------------------------------------#
print "\nterminated !\n";
