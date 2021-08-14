#!/usr/bin/perl
$file=$ARGV[0];
$new=$file;
$new=~s/FPKM/TPM/;
open FILE,"$file";
open NEW,">$new";
$mark=0;
$line=<FILE>;
chomp $line;
print NEW "$line\n";
while($line=<FILE>){
    chomp $line;
    @info=split(/\t/,$line);
    $gene[$mark]=$info[0];
    $sample=scalar(@info);
    for($i=1;$i<$sample;$i++){
        $$i[$mark]=$info[$i];
        $sum[$i]+=$info[$i];
    }
    $mark++;
}
for($i=0;$i<$mark;$i++){
    for($j=1;$j<$sample;$j++){
        if($sum[$j] ne 0){
            $TPM=($$j[$i]/$sum[$j])*1000000;
        }
        else{
            $TPM=0;
        }
            $result.=$TPM."\t";
    }
    chop $result;
    print NEW "$gene[$i]\t$result\n";
    $result='';
}
