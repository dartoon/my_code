#!/usr/bin/perl

$filename="source_list.asc";

open(READFILE,$filename);

$i=0;

while(<READFILE>) {
    chop($_);
    @data=split(" ",$_);
    $name[$i]=$data[0];
    $ra[$i]=$data[1];
    $dec[$i]=$data[2];

    $command="python3 run_HSC_galaxy.py ".$name[$i]." I ".$ra[$i]." ".$dec[$i];

    printf("Running on target %s \n",$name[$i]);
    system($command);
    sleep(1.0);

    $i=$i+1;     
}
close(READFILE);
$N=$i;
